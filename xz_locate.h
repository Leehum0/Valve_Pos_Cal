/**
 * @file    xz_locate.h
 * @brief   磁条安装位置 (X,Z) 识别 —— 多快照联合 + 多起点精搜
 *          目标平台 STM32U575VGT6 (Cortex-M33, 单精度 FPU)
 *
 * 与参考实现 valve_locator.c / xyz_estimation_u5q_3_8_8.m 的一致性硬约束：
 *   - 数据库布局:  idx = (((x-1)*NY+(y-1))*NZ+(z-1))*NS+(s-1)  (全 1 基)
 *   - 平滑:        13 点滑动平均级联 3 遍 (movmean 边界窗口收缩)
 *   - 选峰:        窗口 {60,100}/{140,180}/{220,260}/{300,340} 内
 *                  平滑导数负区主峰 + 三点抛物线细化
 *   - 亚采样插值:  4 点三次卷积 (cubicConv4 等价核)
 *   - 三维插值:    interpXYZ4D_2 语义 (含边界钳位修复)
 *
 * 相对参考实现 vl_calibrate 内部坐标估计段的改进 (依据闭环验证):
 *   1. 多快照联合: 阀门停 K 个行程位置各采一组, X/Z 共享, y_k 独立;
 *   2. 多起点精搜: 粗搜 Top-6 个彼此距离>=2 格的候选盆地分别精搜,
 *      取全局 J 最小 —— 解决曲线族近共线导致的单起点落错盆地问题;
 *   3. 现场偏置提取 xz_extract_offset(): 全行程扫描曲线上用与建库
 *      同源的零位机制取各传感器偏置 (数据库常数偏置误差达 ±0.06,
 *      而 X 识别容忍度约 ±0.01, 必须现场测).
 *
 * 使用流程 (标定期一次性):
 *   1) 全行程慢速扫描 -> xz_extract_offset() x4 -> b_s
 *   2) 阀门停 5 个开度, 每处 ADC 均值 -> xz_normalize() -> 快照
 *   3) xz_estimate() -> X̂/Ẑ + x0/z0/alpha (喂给 vl_build_comp_sig)
 *
 * 无动态内存、无递归。C99。
 */

#ifndef XZ_LOCATE_H
#define XZ_LOCATE_H

#ifdef __cplusplus
extern "C" {
#endif

/* ===================== 数据库维度 (与 valve_locator.h 相同) ========== */
#define XZ_NX            4      /**< 横向采样点数                      */
#define XZ_NY            401    /**< 行程采样点数                      */
#define XZ_NZ            13     /**< 纵向采样点数                      */
#define XZ_NS            4      /**< HMC1501 传感器个数                */

/* ===================== 物理常数 ====================================== */
#define XZ_TRAV_LENGTH   130.12f                          /**< 总行程 mm */
#define XZ_SAMPLE_STEP   (XZ_TRAV_LENGTH / (XZ_NY - 1))   /**< mm/采样点 */
#define XZ_STEP_LATER    1.5f    /**< 横向网格间距 mm                   */
#define XZ_STEP_LONG     0.75f   /**< 纵向网格间距 mm                   */
#define XZ_LONG_OFFSET   8.0f    /**< 纵向起始距离 mm                   */

/* ===================== 算法参数 ====================================== */
#define XZ_MAX_SNAP      8      /**< 快照数上限 (推荐 K=5)             */
#define XZ_N_STARTS      6      /**< 多起点精搜的起点个数              */
#define XZ_START_SEP     2      /**< 起点间最小切比雪夫距离 (格)       */
#define XZ_SMOOTH_LEN    13     /**< 滑动平均窗口 (一致性硬约束)       */
#define XZ_SMOOTH_PASSES 3      /**< 级联遍数     (一致性硬约束)       */
#define XZ_SWEEP_MAX     2048   /**< 偏置提取扫描曲线最大点数          */

/* 行程贴边判定 (主脚本 287 行同源): y<=8 或 y>=8+120/step ≈ 377      */
#define XZ_EDGE_LO       8
#define XZ_EDGE_HI       377

/* 置信度: 次优盆地 J / 最优盆地 J 低于该比值时报 XZ_WARN_AMBIG        */
#define XZ_AMBIG_RATIO   1.5f

/* 快照行程分布下限 (采样点): K 个快照的 y 极差低于该值时报
 * XZ_WARN_SPREAD —— 快照挤在同一开度会退化成单快照可靠性             */
#define XZ_MIN_SPREAD    40.0f

/* ===================== 返回码 ======================================== */
typedef enum {
    XZ_OK           =  0,   /**< 成功                                   */
    XZ_WARN_AMBIG   =  1,   /**< 成功但次优盆地接近, 置信度低           */
    XZ_WARN_SPREAD  =  2,   /**< 成功但快照开度间隔不足, 结果不可靠;
                                 检查阀门是否真的移动了 (bsp_valve_goto) */
    XZ_ERR_PARAM    = -1,   /**< 空指针 / K 越界 / 点数非法 / 信号非有限 */
    XZ_ERR_EDGE     = -2,   /**< 有快照行程序号贴边, 该开度需更换       */
    XZ_ERR_PEAK     = -3    /**< 偏置提取失败 (正峰或峰在窗口边界)      */
} xz_status_t;

/* ===================== 归一化参数 ==================================== */
/**
 * norm[s] = (adc[s] - offset[s]) / scale[s]
 * scale : 建库导出的每传感器半峰峰值 (见 XZ_SCALE_DEFAULT);
 *         数据库重建后必须同步更新。
 * offset: 现场提取 (xz_extract_offset), 勿用出厂常数 —— 52 条建库曲线
 *         的偏置分布跨度达 456~805 ADC 码, 远超 X 识别容忍度。
 */
typedef struct {
    float scale[XZ_NS];
    float offset[XZ_NS];
} xz_norm_t;

/** 20260706 数据库反演的缩放参考值 (参考块 X1Z3, 拟合残差 4e-6) */
#define XZ_SCALE_DEFAULT { 7695.8399f, 7517.6405f, 7284.6545f, 7843.2578f }

/* ===================== 估计结果 ====================================== */
typedef struct {
    float x_idx;                 /**< 横向连续格坐标 (1 基, 1..4)       */
    float z_idx;                 /**< 纵向连续格坐标 (1 基, 1..13)      */
    float lateral_mm;            /**< 横向物理坐标 = (x_idx-1)*1.5      */
    float long_mm;               /**< 纵向物理坐标 = 8+(z_idx-1)*0.75   */
    int   x_near;                /**< 最近标定曲线号 X (1..4)           */
    int   z_near;                /**< 最近标定曲线号 Z (1..13)          */
    int   x0, z0;                /**< 双线性合成左下角整数索引 (1 基)   */
    float alpha_x, alpha_z;      /**< 双线性小数权重 [0,1]              */
    float y_idx[XZ_MAX_SNAP];    /**< 各快照行程序号 (诊断用, 1 基)     */
    float j_best;                /**< 最优盆地终点 J (K 快照累加)       */
    float j_second;              /**< 次优盆地终点 J (置信度评估)       */
} xz_result_t;

/* ===================== API =========================================== */

/** 4 路原始 ADC 均值 -> 归一化信号。scale 含 0/负值时返回 XZ_ERR_PARAM */
xz_status_t xz_normalize(const float adc[XZ_NS], const xz_norm_t *np,
                         float out[XZ_NS]);

/* 偏置合理性校验: 提取的零位幅值须落在扫描曲线中点 ± 该比例 x 半峰峰值
 * 之内 (数据库实测零位幅值距中点最大 0.11 半峰峰); 超出即判选峰错误,
 * 返回 XZ_ERR_PEAK, 调用方应重新扫描。用于拦截速度畸变导致的错误选峰 */
#define XZ_OFFSET_MID_TOL   0.20f

/* 偏置提取选峰窗口的放宽余量 (占扫描长度比例): 容忍扫描速度不均匀
 * 造成的特征位置漂移。窗口按 [建库窗口/401 x n] 缩放后两侧各放宽该值。
 * 上限约束: 建库窗口与相邻信号峰/谷 (正导数区) 的间隙仅 ~5% 行程,
 * 放宽过大会把正导数区圈入窗口导致选峰失败 —— 不得超过 0.05。
 * 对扫描的要求: 速度须足够均匀, 特征位置漂移 <= 该余量 (全程速度
 * 单调漂移 < ~15% 即可满足)。 */
#define XZ_WIN_MARGIN_FRAC  0.04f

/**
 * 全行程扫描曲线中提取一路传感器的偏置 (原始 ADC 单位)。
 * 机制与建库归一化同源: 中心差分导数 -> movmean x3 -> 该传感器
 * 选峰窗口内取负区主峰 -> 抛物线细化 -> cubicConv4 取该处幅值。
 * 两项密度自适应 (建库按 401 点定参, 扫描密度不同必须换算):
 *   - 平滑窗宽 w = 奇数化(13 * n / 401), 下限 13;
 *   - 选峰窗口 = 建库窗口({60,100}等)按 n/401 缩放, 两侧各放宽
 *     XZ_WIN_MARGIN_FRAC (扫描须匀速覆盖全行程, 特征位置按比例出现)。
 * @param sweep      该传感器的扫描曲线 (原始 ADC 值, 等间隔采样,
 *                   从全关到全开单向匀速覆盖全行程)
 * @param n          点数, 32 <= n <= XZ_SWEEP_MAX
 * @param sensor_idx 传感器序号 0..3 (决定选峰窗口)
 * @param offset_raw 输出: 零位幅值 (即该传感器偏置 b_s)
 * @note  多路曲线共用内部静态工作区, 非线程安全(标定期单线程调用)。
 */
xz_status_t xz_extract_offset(const float *sweep, int n, int sensor_idx,
                              float *offset_raw);

/**
 * 多快照多起点 (X,Z) 估计。
 * @param db     四维数据库 (Flash), 布局见文件头
 * @param snaps  K 组已归一化的 4 路信号
 * @param k      快照数, 1 <= k <= XZ_MAX_SNAP (推荐 5; k=1 兼容参考实现
 *               的单快照模式, 仅建议行程中段使用)
 * @param res    输出结果
 * @return XZ_OK / XZ_WARN_* (结果已填) / XZ_ERR_*
 * @note  内部使用静态候选缓冲, 非线程安全 (标定期单线程调用)。
 */
xz_status_t xz_estimate(const float *db, const float snaps[][XZ_NS], int k,
                        xz_result_t *res);

#ifdef __cplusplus
}
#endif
#endif /* XZ_LOCATE_H */
