/**
 * @file    xz_locate.c
 * @brief   磁条安装位置 (X,Z) 识别 —— 实现 (见 xz_locate.h)
 *
 * 逐函数对应关系:
 *   xz_db_at1          <- valve_locator.c vl_db_at1 (功能对应, 但布局改为
 *                         与归一化 txt 行序一致的 x->z->y->s; 若与
 *                         valve_locator.c 集成, vl_db_at1 须同步改公式)
 *   xz_cubic_conv4     <- cubicConv4.m
 *   xz_interp          <- interpXYZ4D_2.m (含边界钳位修复)
 *   xz_movmean/deriv   <- movmean(13)x3 (build_database_18_9 同源)
 *   xz_extract_offset  <- 建库归一化的零位幅值机制 (readme.txt)
 *   xz_estimate        <- estimateXZRatio_7.m 扩展:
 *                         多快照联合 + 多起点精搜 (改进依据见 .h 头注释)
 *
 * 算法经 20260706 数据库闭环验证 (Python 同构镜像):
 *   52 网格点自识别 52/52; 亚网格 K=5 无扰动 max|Xerr|=0.27 格,
 *   现实误差模型(偏置±0.005/增益±2%/噪声σ0.01) max|Xerr|≈1.1 格,
 *   max|Zerr|<=0.5 格。详见《磁条安装位置XZ识别方法_v1_0.md》。
 */

#include <stddef.h>
#include <math.h>
#include "xz_locate.h"

/* ================= 展平数据库按 1 基 (x,y,z,s) 取值 ================= */
/* 布局与归一化 txt 行序一致 (x 外层 -> z 中层 -> y 内层 -> s):
 * 沿 y 扫描为连续内存访问, 对 Flash 预取/缓存友好 */
static float xz_db_at1(const float *db, int x, int y, int z, int s)
{
    return db[(((x - 1) * XZ_NZ + (z - 1)) * XZ_NY + (y - 1)) * XZ_NS + (s - 1)];
}

/* ===================== cubicConv4.m ================================= */
static float xz_cubic_conv4(float pm1, float p0, float p1, float p2, float t)
{
    return p0 + 0.5f * t * ((p1 - pm1)
        + t * (2.0f * pm1 - 5.0f * p0 + 4.0f * p1 - p2)
        + t * t * (3.0f * (p0 - p1) + p2 - pm1));
}

/* ============ interpXYZ4D_2: y 三次卷积, x/z 双线性, 1 基 ============ */
static void xz_interp(const float *db, float x, float y, float z,
                      float out[XZ_NS])
{
    int   x0, y0, z0, s;
    float ax, az, ty;

    if (x < 1.0f)               { x = 1.0f; }
    if (x > (float)XZ_NX)       { x = (float)XZ_NX; }
    if (z < 1.0f)               { z = 1.0f; }
    if (z > (float)XZ_NZ)       { z = (float)XZ_NZ; }
    if (y < 2.0f)               { y = 2.0f; }          /* 4 邻点模板 */
    if (y > (float)(XZ_NY - 2)) { y = (float)(XZ_NY - 2); }

    x0 = (int)floorf(x); if (x0 > XZ_NX - 1) { x0 = XZ_NX - 1; } ax = x - (float)x0;
    z0 = (int)floorf(z); if (z0 > XZ_NZ - 1) { z0 = XZ_NZ - 1; } az = z - (float)z0;
    y0 = (int)floorf(y); ty = y - (float)y0;

    for (s = 1; s <= XZ_NS; s++) {
        float c[4];
        int   q;
        static const int DXZ[4][2] = { {0,0}, {1,0}, {0,1}, {1,1} };
        for (q = 0; q < 4; q++) {
            int xi = x0 + DXZ[q][0], zi = z0 + DXZ[q][1];
            c[q] = xz_cubic_conv4(xz_db_at1(db, xi, y0 - 1, zi, s),
                                  xz_db_at1(db, xi, y0,     zi, s),
                                  xz_db_at1(db, xi, y0 + 1, zi, s),
                                  xz_db_at1(db, xi, y0 + 2, zi, s), ty);
        }
        out[s - 1] = c[0] * (1 - ax) * (1 - az) + c[1] * ax * (1 - az)
                   + c[2] * (1 - ax) * az       + c[3] * ax * az;
    }
}

/* ===== 居中滑动平均, 边界窗口收缩 (等价 MATLAB movmean), 窗宽 w ===== */
static void xz_movmean(const float *x, float *y, int n, int w)
{
    const int half = w / 2;
    int i;
    for (i = 0; i < n; i++) {
        int   a = i - half, b = i + half, j;
        float acc = 0.0f;
        if (a < 0)     { a = 0; }
        if (b > n - 1) { b = n - 1; }
        for (j = a; j <= b; j++) { acc += x[j]; }
        y[i] = acc / (float)(b - a + 1);
    }
}

/* ================== 归一化 ========================================== */
xz_status_t xz_normalize(const float adc[XZ_NS], const xz_norm_t *np,
                         float out[XZ_NS])
{
    int s;
    if (adc == NULL || np == NULL || out == NULL) { return XZ_ERR_PARAM; }
    for (s = 0; s < XZ_NS; s++) {
        if (!(np->scale[s] > 0.0f)) { return XZ_ERR_PARAM; }  /* 含 NaN */
        out[s] = (adc[s] - np->offset[s]) / np->scale[s];
    }
    return XZ_OK;
}

/* ================== 现场偏置提取 ==================================== */
/* 标定期单线程调用, 静态工作区 2 x 8 KB */
static float xz_wk1[XZ_SWEEP_MAX];
static float xz_wk2[XZ_SWEEP_MAX];

/* 各传感器建库选峰窗口 (1 基, 401 点坐标系), 与 build_database_18_9 一致 */
static const int XZ_WIN[XZ_NS][2] = {
    {  60, 100 }, { 140, 180 }, { 220, 260 }, { 300, 340 }
};

xz_status_t xz_extract_offset(const float *sweep, int n, int sensor_idx,
                              float *offset_raw)
{
    int   i, p, pk, w, lo, hi, margin;
    float pv, delta, denom;

    if (sweep == NULL || offset_raw == NULL || n < 32 || n > XZ_SWEEP_MAX
        || sensor_idx < 0 || sensor_idx >= XZ_NS) {
        return XZ_ERR_PARAM;
    }

    /* 平滑窗口随采样密度自适应: 与建库 (401 点 / 13 点窗) 覆盖相同的
     * 物理跨度。n=401 -> 13, n=1200 -> 39。奇数化, 下限 13。 */
    w = (XZ_SMOOTH_LEN * n) / XZ_NY;
    if (w < XZ_SMOOTH_LEN) { w = XZ_SMOOTH_LEN; }
    w |= 1;

    /* 中心差分导数 (gradient 同源) */
    xz_wk1[0]     = sweep[1] - sweep[0];
    xz_wk1[n - 1] = sweep[n - 1] - sweep[n - 2];
    for (i = 1; i < n - 1; i++) {
        xz_wk1[i] = 0.5f * (sweep[i + 1] - sweep[i - 1]);
    }
    /* movmean(w) 级联 3 遍, 乒乓缓冲 */
    for (p = 0; p < XZ_SMOOTH_PASSES; p++) {
        if (p & 1) { xz_movmean(xz_wk2, xz_wk1, n, w); }
        else       { xz_movmean(xz_wk1, xz_wk2, n, w); }
    }
    /* 选峰窗口: 建库窗口按扫描密度缩放 + 两侧放宽 (0 基) */
    margin = (int)(XZ_WIN_MARGIN_FRAC * (float)n);
    lo = (XZ_WIN[sensor_idx][0] * n) / XZ_NY - margin;
    hi = (XZ_WIN[sensor_idx][1] * n) / XZ_NY + margin;
    if (lo < 1)     { lo = 1; }
    if (hi > n - 2) { hi = n - 2; }
    if (hi - lo < 4) { return XZ_ERR_PARAM; }

    {
        const float *d = (XZ_SMOOTH_PASSES & 1) ? xz_wk2 : xz_wk1;

        /* 窗口内负区主峰 (导数最大值, 须为负且不在窗口边界):
         * 与 build_database_18_9 的选峰语义同源 */
        pk = lo; pv = d[lo];
        for (i = lo + 1; i <= hi; i++) {
            if (d[i] > pv) { pv = d[i]; pk = i; }
        }
        if (pv >= 0.0f || pk == lo || pk == hi) { return XZ_ERR_PEAK; }

        /* 三点抛物线细化 */
        denom = d[pk - 1] - 2.0f * d[pk] + d[pk + 1];
        delta = (denom != 0.0f) ? 0.5f * (d[pk - 1] - d[pk + 1]) / denom : 0.0f;
    }

    /* 4 点三次卷积取亚采样位置幅值 (需 pk-1..pk+2 有效) */
    if (delta < 0.0f) { pk -= 1; delta += 1.0f; }
    if (pk < 1)       { pk = 1; delta = 0.0f; }
    if (pk > n - 3)   { pk = n - 3; delta = 1.0f; }
    {
        float off = xz_cubic_conv4(sweep[pk - 1], sweep[pk],
                                   sweep[pk + 1], sweep[pk + 2], delta);
        /* 合理性校验: 零位幅值应在曲线中点附近, 拦截速度畸变下的错误选峰 */
        float mx = sweep[0], mn = sweep[0];
        for (i = 1; i < n; i++) {
            if (sweep[i] > mx) { mx = sweep[i]; }
            if (sweep[i] < mn) { mn = sweep[i]; }
        }
        if (fabsf(off - 0.5f * (mx + mn))
                > XZ_OFFSET_MID_TOL * 0.5f * (mx - mn)) {
            return XZ_ERR_PEAK;
        }
        *offset_raw = off;
    }
    return XZ_OK;
}

/* ================== (X,Z) 估计 ====================================== */

/* 粗搜候选: 每个 (x,z) 网格的 K 快照联合 J 与各快照最优 y */
typedef struct {
    float j;
    int   x, z;
    int   y[XZ_MAX_SNAP];
} xz_cand_t;

/* 单快照在整数网格 (x,z) 的全行程扫描: 返回最小 J, *y_best 收最优行程 */
static float xz_scan_y(const float *db, const float sen[XZ_NS],
                       int x, int z, int *y_best)
{
    int   y, s, yb = 1;
    float jb = 1e30f;
    for (y = 1; y <= XZ_NY; y++) {
        float j = 0.0f;
        for (s = 1; s <= XZ_NS; s++) {
            j += fabsf(xz_db_at1(db, x, y, z, s) - sen[s - 1]);
            if (j >= jb) { break; }        /* 已劣于当前最优, 提前放弃 */
        }
        if (j < jb) { jb = j; yb = y; }
    }
    *y_best = yb;
    return jb;
}

/* 候选 (x,z) 处的联合 J: 每快照在 y_init[k]±2h 内 5 点取最小, y_out 收各自最优 */
static float xz_eval(const float *db, const float snaps[][XZ_NS], int k,
                     float x, float z, const float *y_init, float h,
                     float *y_out)
{
    int   i, dy, s;
    float jtot = 0.0f;
    for (i = 0; i < k; i++) {
        float jb = 1e30f, yb = y_init[i];
        for (dy = -2; dy <= 2; dy++) {
            float y = y_init[i] + (float)dy * h, v[XZ_NS], j = 0.0f;
            if (y < 1.0f || y > (float)XZ_NY) { continue; }
            xz_interp(db, x, y, z, v);
            for (s = 0; s < XZ_NS; s++) { j += fabsf(v[s] - snaps[i][s]); }
            if (j < jb) { jb = j; yb = y; }
        }
        jtot += jb;
        y_out[i] = yb;
    }
    return jtot;
}

/* 单起点三层由粗到细精搜 (estimateXZRatio_7 的多快照扩展) */
static float xz_fine(const float *db, const float snaps[][XZ_NS], int k,
                     const xz_cand_t *start, float *x_out, float *z_out,
                     float *y_out)
{
    static const float steps[3] = { 0.5f, 0.125f, 0.03125f };
    float xb = (float)start->x, zb = (float)start->z;
    float yb[XZ_MAX_SNAP], jb = 1e30f;
    int   lv, dx, dz, i;

    for (i = 0; i < k; i++) { yb[i] = (float)start->y[i]; }

    for (lv = 0; lv < 3; lv++) {
        float h = steps[lv];
        float jbest = 1e30f, x1 = xb, z1 = zb;
        float y1[XZ_MAX_SNAP], yc[XZ_MAX_SNAP];
        for (i = 0; i < k; i++) { y1[i] = yb[i]; }

        for (dx = -2; dx <= 2; dx++) {
            float x = xb + (float)dx * h;
            if (x < 1.0f || x > (float)XZ_NX) { continue; }
            for (dz = -2; dz <= 2; dz++) {
                float z = zb + (float)dz * h, j;
                if (z < 1.0f || z > (float)XZ_NZ) { continue; }
                j = xz_eval(db, snaps, k, x, z, yb, h, yc);
                if (j < jbest) {
                    jbest = j; x1 = x; z1 = z;
                    for (i = 0; i < k; i++) { y1[i] = yc[i]; }
                }
            }
        }
        xb = x1; zb = z1; jb = jbest;
        for (i = 0; i < k; i++) { yb[i] = y1[i]; }
    }

    *x_out = xb; *z_out = zb;
    for (i = 0; i < k; i++) { y_out[i] = yb[i]; }
    return jb;
}

xz_status_t xz_estimate(const float *db, const float snaps[][XZ_NS], int k,
                        xz_result_t *res)
{
    /* 静态候选缓冲 (~2.3 KB): 避免压 RTOS 任务栈; 标定期单线程调用 */
    static xz_cand_t cand[XZ_NX * XZ_NZ];
    static int       order[XZ_NX * XZ_NZ];
    int       starts[XZ_N_STARTS];
    int       n_starts = 0;
    int       x, z, i, j, c;
    float     best_j = 1e30f, best_x = 1.0f, best_z = 1.0f;
    float     best_y[XZ_MAX_SNAP];
    float     second_j = 1e30f;

    if (db == NULL || snaps == NULL || res == NULL
        || k < 1 || k > XZ_MAX_SNAP) {
        return XZ_ERR_PARAM;
    }
    /* 信号有限性校验: NaN/inf 会让比较全部失效并静默给出 (1,1) */
    for (i = 0; i < k; i++) {
        for (j = 0; j < XZ_NS; j++) {
            if (!isfinite(snaps[i][j])) { return XZ_ERR_PARAM; }
        }
    }

    /* ---- 1. 粗搜: 全部 4x13 网格, 每快照全行程扫描 ---- */
    c = 0;
    for (x = 1; x <= XZ_NX; x++) {
        for (z = 1; z <= XZ_NZ; z++) {
            cand[c].x = x; cand[c].z = z; cand[c].j = 0.0f;
            for (i = 0; i < k; i++) {
                cand[c].j += xz_scan_y(db, snaps[i], x, z, &cand[c].y[i]);
            }
            order[c] = c;
            c++;
        }
    }

    /* 按 J 升序 (插入排序, 52 项) */
    for (i = 1; i < c; i++) {
        int t = order[i];
        for (j = i - 1; j >= 0 && cand[order[j]].j > cand[t].j; j--) {
            order[j + 1] = order[j];
        }
        order[j + 1] = t;
    }

    /* ---- 2. 选 XZ_N_STARTS 个彼此距离 >= XZ_START_SEP 格的起点 ---- */
    for (i = 0; i < c && n_starts < XZ_N_STARTS; i++) {
        const xz_cand_t *ci = &cand[order[i]];
        int ok = 1;
        for (j = 0; j < n_starts; j++) {
            const xz_cand_t *cj = &cand[starts[j]];
            int ddx = ci->x - cj->x, ddz = ci->z - cj->z;
            if (ddx < 0) { ddx = -ddx; }
            if (ddz < 0) { ddz = -ddz; }
            if ((ddx > ddz ? ddx : ddz) < XZ_START_SEP) { ok = 0; break; }
        }
        if (ok) { starts[n_starts++] = order[i]; }
    }

    /* ---- 3. 每起点三层精搜, 取全局 J 最小; 记录竞争盆地 J ---- */
    for (i = 0; i < n_starts; i++) {
        float xf, zf, yf[XZ_MAX_SNAP];
        float jf = xz_fine(db, snaps, k, &cand[starts[i]], &xf, &zf, yf);
        if (jf < best_j) {
            /* 原最优若与新最优不同盆地, 降级为次优 */
            {
                float ddx = fabsf(best_x - xf), ddz = fabsf(best_z - zf);
                if (best_j < 1e29f && (ddx >= 1.0f || ddz >= 1.0f)
                    && best_j < second_j) {
                    second_j = best_j;
                }
            }
            best_j = jf; best_x = xf; best_z = zf;
            for (j = 0; j < k; j++) { best_y[j] = yf[j]; }
        } else {
            float ddx = fabsf(best_x - xf), ddz = fabsf(best_z - zf);
            if ((ddx >= 1.0f || ddz >= 1.0f) && jf < second_j) {
                second_j = jf;
            }
        }
    }

    /* ---- 4. 填充结果 ---- */
    if (best_x < 1.0f)             { best_x = 1.0f; }
    if (best_x > (float)XZ_NX)     { best_x = (float)XZ_NX; }
    if (best_z < 1.0f)             { best_z = 1.0f; }
    if (best_z > (float)XZ_NZ)     { best_z = (float)XZ_NZ; }

    res->x_idx      = best_x;
    res->z_idx      = best_z;
    res->lateral_mm = (best_x - 1.0f) * XZ_STEP_LATER;
    res->long_mm    = XZ_LONG_OFFSET + (best_z - 1.0f) * XZ_STEP_LONG;
    res->x_near     = (int)floorf(best_x + 0.5f);
    res->z_near     = (int)floorf(best_z + 0.5f);
    res->j_best     = best_j;
    res->j_second   = second_j;

    /* x0/z0/alpha (estimateXZRatio_7 输出约定, 供 vl_build_comp_sig) */
    res->x0 = (int)floorf(best_x);
    if (res->x0 >= XZ_NX) { res->x0 = XZ_NX - 1; res->alpha_x = 1.0f; }
    else                  { res->alpha_x = best_x - (float)res->x0; }
    res->z0 = (int)floorf(best_z);
    if (res->z0 >= XZ_NZ) { res->z0 = XZ_NZ - 1; res->alpha_z = 1.0f; }
    else                  { res->alpha_z = best_z - (float)res->z0; }
    if (res->alpha_x < 0.0f) { res->alpha_x = 0.0f; }
    if (res->alpha_x > 1.0f) { res->alpha_x = 1.0f; }
    if (res->alpha_z < 0.0f) { res->alpha_z = 0.0f; }
    if (res->alpha_z > 1.0f) { res->alpha_z = 1.0f; }

    for (i = 0; i < k; i++)          { res->y_idx[i] = best_y[i]; }
    for (i = k; i < XZ_MAX_SNAP; i++) { res->y_idx[i] = 0.0f; }

    /* ---- 5. 有效性 ---- */
    for (i = 0; i < k; i++) {
        if (best_y[i] <= (float)XZ_EDGE_LO || best_y[i] >= (float)XZ_EDGE_HI) {
            return XZ_ERR_EDGE;
        }
    }
    /* 快照行程分布: 挤在同一开度会退化成单快照可靠性 */
    if (k >= 2) {
        float ymin = best_y[0], ymax = best_y[0];
        for (i = 1; i < k; i++) {
            if (best_y[i] < ymin) { ymin = best_y[i]; }
            if (best_y[i] > ymax) { ymax = best_y[i]; }
        }
        if (ymax - ymin < XZ_MIN_SPREAD) {
            return XZ_WARN_SPREAD;
        }
    }
    if (second_j < best_j * XZ_AMBIG_RATIO) {
        return XZ_WARN_AMBIG;
    }
    return XZ_OK;
}
