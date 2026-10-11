# SDM 弯曲断层棋盘测试说明

本文说明 [make_sdm_checkerboard.py](make_sdm_checkerboard.py) 的网格分配，以及 [run_sdm_checkerboard.py](run_sdm_checkerboard.py) 的自动化正演、反演和恢复指标计算。脚本仅依赖 NumPy，使用已有 SDM2025 执行文件，不修改科学内核。

## 1. Maduo 测试设置

文中六组 GPS、150 条观测及末尾数值记录的是 2026-10-10 首次验证。当前模板已增加两组 InSAR，实际数据集由 `ngd` 指定，程序始终处理全部数据集。

程序读取 [sdm2025_template.inp](../Maduo/Checkboard/sdm2025_template.inp) 和 [checkerboard_15km.dat](../Maduo/Checkboard/checkerboard_15km.dat)。模板指定的原反演模型 `Maduo/Checkboard/slip_model.dat` 用于核对几何、patch 数目和行序；已复制的 `ZHYRELAX.txt` 原样用于反演。

**现有棋盘实际为走向 15 km × 倾向 9 km、滑动 3 m / 0 m、非零区 rake 15°。文件名中的 15 km 不表示两个方向均为 15 km。自动化程序使用这个模型，不重新生成棋盘。**

| 断层段 | patch 数目 | 走向 × 倾向棋盘分箱数 | 非零滑动 patch | 零滑动 patch |
| --- | ---: | --- | ---: | ---: |
| 1 | 741 | 12 × 5 | 372 | 369 |
| 2 | 169 | 3 × 5 | 86 | 83 |
| 合计 | 910 | 每段独立分配 | 458 | 452 |

本模板 `iearth=0`、泊松比 0.25，使用均匀弹性半空间中的 Okada 矩形位错解，线性叠加各 patch 的静态位移。两段 rake 范围均为 −30°～30°，最大滑动均为 7 m。测试保留地球模型、断层几何、观测分布、单位、方向、误差和反演约束，评估这些假设下的滑动恢复能力。

## 2. 自动运行全部数据集

在仓库根目录运行无噪声测试：

~~~bash
python Script/run_sdm_checkerboard.py \
    Maduo/Checkboard/sdm2025_template.inp \
    Maduo/Checkboard/checkerboard_15km.dat
~~~

加入一倍原观测误差的高斯噪声：

~~~bash
python Script/run_sdm_checkerboard.py \
    Maduo/Checkboard/sdm2025_template.inp \
    Maduo/Checkboard/checkerboard_15km.dat \
    --noise-scale 1 --seed 0
~~~

| 参数 | 含义 | 默认值 |
| --- | --- | --- |
| 第一个位置参数 | 原反演模板，须有正迭代次数和已有滑动输出模型 | 必选 |
| 第二个位置参数 | 棋盘真值，须有一行表头和明确的分段空行 | 必选 |
| `--sdm PATH` | SDM2025 执行文件 | 仓库内 `SDM2025/SourceCode/sdm2025` |
| `--noise-scale K` | 有限且非负的高斯噪声倍率 | `0` |
| `--seed N` | 非负整数随机种子 | `0` |

模板内部的相对路径均从 **模板所在目录** 解析，SDM 也在该目录运行；所有生成文件写入该目录。启动 Python 时所在目录只影响命令行参数中的相对路径。本例 `../GPS/...` 指向 `Maduo/GPS/`。

程序处理 `ngd` 指定的全部数据集：

| 顺序 | 原观测文件 basename | 观测数 | 模板权重 |
| --- | --- | ---: | ---: |
| 1 | `coseis_maduo2021all.ns1` | 44 | 1 |
| 2 | `coseis_maduo2021all.ew1` | 44 | 1 |
| 3 | `coseis_maduo2021all.ns2` | 4 | 2 |
| 4 | `coseis_maduo2021all.ew2` | 4 | 2 |
| 5 | `coseis_maduotri.ns` | 27 | 5 |
| 6 | `coseis_maduotri.ew` | 27 | 5 |
| 合计 | 六组全部参与正演与反演 | 150 | SDM 中归一化 |

固定方向文件为四列 `lat lon observation sigma`；逐点方向文件为六列 `lat lon observation sigma incidence azimuth`。`nheader` 指定表头物理行数。方向标志为 1 时使用模板的固定角度，否则使用文件中的逐点角度。数据区每个测点须占一个非空行，误差严格大于零，所有数值有限。

## 3. 配置保护、命名与重跑

程序依照 [sdmgetinp.f](../SDM2025/SourceCode/sdmgetinp.f) 的 Fortran 记录顺序提取非空、非注释记录。保留地球模型、断层几何、离散化、滑动约束、`ngd/datunit/nheader`、数据集权重与方向、平滑设置和 Zhang 开关。生成配置不含空行；自行添加注释时须以 `#` 开头。

只调整 correction、迭代次数及流程所需文件名：

1. 开始前检查执行文件、观测、已有反演模型、地球模型文件、patch 文件（`idisc=0` 时）、松弛因子和所有生成文件名。逐段核对棋盘与已有反演模型的几何、patch 数目及行序。文件名碰撞、覆盖原文件、输出符号链接或硬链接立即报错。
2. 清理本次测试的旧生成文件，重新复制原始棋盘。不能复用上次已被正演重写的工作模型。
3. 正演迭代次数设为零，所有数据集生成预测；SDM 覆盖工作模型并重算应力列。
4. 同时检查退出码、正常完成标志、模型和全部输出文件。Fortran `STOP` 可能返回零，不能仅依赖退出码。检查几何、滑动、观测数及坐标行序；坐标允许五位小数输出造成的约 `5×10⁻⁶` 度舍入误差。
5. 逐文件读取正演第 5 列 prediction，加入可选噪声，另存合成观测。**仅替换原观测第 3 列数字文本**；表头、坐标、误差、逐点方向、空白字符、换行和行序保留。拟合输出第 4 列是残差，不能当成误差。
6. 反演引用全部合成观测，恢复模板原迭代次数，本例为 10000。保留原松弛文件路径、平滑 `1 0.005 0` 及其它科学设置。检查恢复模型及全部预测后计算指标。

| 用途 | 本例文件名 |
| --- | --- |
| 正演配置 | `work_sdm2025_template.inp` |
| 反演配置 | `invwork_sdm2025_template.inp` |
| 工作棋盘 / 正演模型 | `work_checkerboard_15km.dat` |
| 恢复模型 | `invwork_checkerboard_15km.dat` |
| 六组正演预测 | `work_far_ns1.dat`、`work_far_ew1.dat`、`work_far_ns2.dat`、`work_far_ew2.dat`、`work_tri_ns.dat`、`work_tri_ew.dat` |
| 六组恢复预测 | 在上述文件名前直接加 `inv`，如 `invwork_far_ns1.dat` |
| 六组合成观测 | 原观测 basename 前加 `input_`，如 `input_coseis_maduo2021all.ns1` |
| 两阶段控制台记录 | `work_forward.log`、`invwork_inverse.log` |
| 模型恢复指标 | `invwork_resolution_metrics.csv` |
| 位移拟合指标 | `invwork_dataset_fit.csv` |

SDM 还会写 `s01_work_...`、`s02_work_...`、`s01_invwork_...`、`s02_invwork_...` 分段模型和 `log_work_...`、`log_invwork_...` 原生迭代日志。

**重复运行自动覆盖上述本次测试文件。** 如需比较不同噪声或棋盘试验，应先另存本次生成文件。任何阶段失败立即停止；上次的合成观测、恢复模型和 CSV 已在开始时清理，不能误用。失败阶段的配置和控制台记录保留供检查。原始模板、棋盘、观测、松弛文件及 `Maduo/Inv/` 结果只读。

## 4. correction 参数如何设置

[sdmiterate.f](../SDM2025/SourceCode/sdmiterate.f) 在 `niter=0` 时将全部 `corrpar` 清零；随后 [sdmdatfit.f](../SDM2025/SourceCode/sdmdatfit.f) 计算：

~~~text
断层位移：datmdl  = G * s
数据修正：corrmdl = C * p
总拟合量：G * s + C * p
~~~

`s` 是内部走滑、上倾滑动向量，`C` 是修正矩阵，`p` 是 `corrpar`。参数上下限不是正演参数值；将上下限设成相同的非零数，也不能指定正演修正。程序不会读入已有 `corpara.dat` 作为正演参数。

**本自动化测试正演、反演均设 `ndpar=0`，不考虑修正项。** 程序将模板中的：

~~~text
5 './offset_euler_gps.grn'
-0.5 0.5
-0.5 0.5
-5.0E-07 5.0E-07
-5.0E-07 5.0E-07
-5.0E-07 5.0E-07
~~~

整体替换成一行 `0`，删除五组上下限以及输出列表中的 `corpara.dat`。原模板不变。这张修正矩阵不再是依赖，即使 Checkboard 目录没有它也可运行。

仅把上下限设成 `0 0` 不等于关闭 correction，`ndpar=5` 时 SDM 仍会读取并归一化矩阵。保留 correction 的正演仍输出零修正，而且日志目标函数可能受初始化 correction 分量影响，因此本流程直接关闭这些未知量。

拟合输出六列为纬度、经度、观测、残差、prediction、correction。本流程两阶段第六列均为零，合成观测采用第五列。另行设计非零修正测试时需显式外部加入 `C*p_true`，再另行配置联合反演，不属于本程序。修正矩阵的位移单位为米，加回文件单位时应除以 `datunit`。

[make_offset_euler.py](make_offset_euler.py) 当前生成北、东、上三组的六参数矩阵，不能直接替代本模板六组水平观测所需的五参数矩阵。使用前须明确单位、列含义及数据集和测点行序。

## 5. 高斯噪声与复现

原误差为 `sigma_i`，定义：

~~~text
epsilon_i ~ N(0, (K * sigma_i)^2)
synthetic_i = prediction_i + epsilon_i
~~~

`K` 为 `--noise-scale`。噪声和预测使用原文件单位，不额外转换。Maduo `datunit=0.01`，文件单位为 cm；`sigma_i=0.2`、`K=1` 即噪声标准差 0.2 cm。原误差列不变，反演权重和其它科学设置不变。

使用一个 `numpy.random.default_rng(seed)`，按模板数据集顺序、文件内测点顺序生成独立噪声。`K=0` 直接使用 prediction，不消耗随机数。相同输入、倍率和种子生成相同观测；更改种子得到另一试验。两个 CSV 均记录 `noise_scale` 和 `seed`。这里假设测点噪声独立，不模拟空间相关噪声。

## 6. 面积加权模型恢复指标

以原始棋盘为真值，逐 patch 对齐恢复模型。`A_i=length_km×width_km` 为真值面积，单位 km²。滑动向量使用第 8、9 列走滑、下倾分量，单位 m；幅值重新计算 `hypot(strike,ddip)`，不使用第 10 列。

真值向量为 `s_i`，恢复为 `r_i`；真值、恢复幅值为 `S_i`、`R_i`。非零区 `H={i:S_i>0}`，零区 `L={i:S_i=0}`，不设人为阈值。

| CSV 字段 | 定义 | 理想值 |
| --- | --- | --- |
| `vector_rmse_m` | `sqrt(Σ(A_i·‖r_i−s_i‖²) / ΣA_i)` | 0 m |
| `relative_vector_error` | `sqrt(Σ(A_i·‖r_i−s_i‖²) / Σ(A_i·‖s_i‖²))` | 0 |
| `amplitude_correlation` | 真值与恢复幅值的面积加权 Pearson 相关 | 1 |
| `nonzero_amplitude_recovery` | `Σ_H(A_i·R_i) / Σ_H(A_i·S_i)` | 1 |
| `zero_leakage_rms_m` | `sqrt(Σ_L(A_i·R_i²) / Σ_L A_i)` | 0 m |
| `normalized_leakage` | 泄漏 RMS 除以非零真值区面积加权平均幅值 | 0 |

相关系数为：

~~~text
S_bar = Σ(A*S)/ΣA
R_bar = Σ(A*R)/ΣA
correlation = Σ[A*(S-S_bar)*(R-R_bar)]
              / sqrt(Σ[A*(S-S_bar)^2] * Σ[A*(R-R_bar)^2])
~~~

`invwork_resolution_metrics.csv` 输出全断层 `scope=all` 和各段 `scope=segment`，含段序号、patch 数目、面积、噪声倍率和种子。分母为零、缺少相应区域或幅值方差为零时，英文控制台标记 `N/A`（不可计算），CSV 留空。幅值指标无法识别方向错误，应与向量误差一起解读。

## 7. 位移拟合 WRMS

残差为合成观测减恢复预测，乘 `datunit` 转成米。对数据集 `j`：

~~~text
q_i = sigma_i^(-2)
WRMS_j = sqrt(Σ_j(q_i * residual_m_i^2) / Σ_j q_i)
alpha_j = dataset_weight_j / Σ dataset_weight
WRMS_all = sqrt(Σ_j(alpha_j * WRMS_j^2))
~~~

这与 [sdmdata.f](../SDM2025/SourceCode/sdmdata.f) 一致：组内归一化逆方差，再按归一化数据集权重联合。不能直接混合所有测点逆方差；`datunit` 对同组权重比例的影响相消。

`invwork_dataset_fit.csv` 包含全部数据集的名称、顺序、观测数、归一化权重、WRMS 及联合统计 `scope=all`，位移单位统一为 m。**位移拟合良好不能单独证明滑动分辨率高。** 还应检查幅值恢复、零区泄漏、向量误差、深浅部及边缘结构。几何和地球模型误差需另外设计测试。

## 8. 另行生成棋盘

以下命令另存一个 9 km × 9 km、1 m / 0 m、rake 0° 模型，供其它测试使用；自动化程序不执行此命令：

~~~bash
python Script/make_sdm_checkerboard.py \
    Maduo/Inv/slip_model.dat Maduo/Checkboard/checkerboard_9km.dat \
    --cell-size 9 9 --high-slip 1 --low-slip 0 --rake 0
~~~

| 生成器参数 | 含义 | 默认值 |
| --- | --- | --- |
| 输入、输出位置参数 | 至少 13 列 SDM 模型，输出须与输入不同 | 必选 |
| `--skiprows N` | 跳过并保留前 N 个物理行，不解释表头文字 | `1` |
| `--cell-size DX DY` | 走向、倾向边长，单位 km | `9 9` |
| `--high-slip`、`--low-slip` | 幅值，单位 m，须 `0<=low<high` | `1`、`0` |
| `--rake` | 非零区共同 rake，单位度 | `0` |
| `--origin X0 Y0` | 每段的局部棋盘原点，单位 km | `0 0` |
| `--phase` | 0 为偶数单元高滑动，1 交换高低 | `0` |

模型前 13 列：

~~~text
lat_deg lon_deg depth_km x_local_km y_local_km length_km width_km
slp_strk_m slp_ddip_m slp_am_m strike_deg dip_deg rake_deg
~~~

支持额外数值列，所有行列数须一致。数据区空行分段。无表头用 `--skiprows 0`，多行表头指定实际数目；**直接用于 SDM 或自动化流程必须恰好一行表头**。生成器不读取输入约束，应确认幅值和 rake 满足反演先验。

[sdmdisc.f](../SDM2025/SourceCode/sdmdisc.f) 中 `x_local_km=pl/KM2M` 是 patch 中心沿弯曲顶边迹线的累计距离，`y_local_km=pw/KM2M` 是沿倾向距离，不是垂直深度。每段独立分配：

~~~text
i = floor((x_local_km-X0)/DX)
j = floor((y_local_km-Y0)/DY)
~~~

默认 `i+j` 偶数赋高滑动，奇数赋低滑动，沿两个方向交替。按 patch 中心、左闭右开边界归属；边缘可截短，不移动或拆分 patch。默认原点 `0 0` 为自动离散断层的左上边界。目标边长不强制等于三个 patch 的精确尺寸；原网格约 3 km × 3 km，只能近似表达棋盘边界。更小尺度需先另行离散化。

文件符号约定来自 [sdmoutput.f](../SDM2025/SourceCode/sdmoutput.f)：

~~~text
slp_strk_m = S*cos(rake)
slp_ddip_m = -S*sin(rake)
~~~

内部第二分量为上倾滑动，输出第二分量为下倾滑动，故带负号。生成器只更新这两列、幅值及 rake，保留八位小数；其它文本、几何、行序和分段不变。零幅值 rake 标为 0。**保留的应力列仍来自原模型，需通过正演重算。**

## 9. 验证与无噪声结果

~~~bash
python -m unittest discover -s Script -p test_run_sdm_checkerboard.py
~~~

测试覆盖全部数据集、配置保护、四/六列和表头换行保护、舍入和错位、噪声标准差及种子复现、面积权重手算、方向错误、零区泄漏、不可计算项、文件冲突、原文件保护、零退出却未完成、缺失输出和失败重跑旧结果清理。小型原生 SDM 样例包含四个 patch、两个数据集、四条观测，完成带噪声正演与反演；缺少执行文件时跳过原生样例。

2026-10-10 对当前 Maduo 完成无噪声运行（`K=0, seed=0`），六组 150 条观测全部参与，模板、棋盘、观测、松弛文件及 Inv 文件等 43 个原始文件的 SHA-256 未改变。

| 全断层指标 | 数值 |
| --- | ---: |
| 滑动向量 RMSE | 1.50257 m |
| 滑动向量相对误差 | 0.705985 |
| 幅值相关系数 | 0.233303 |
| 非零区幅值恢复率 | 0.545739 |
| 零区泄漏 RMS | 1.48088 m |
| 归一化泄漏 | 0.493627 |
| 联合位移 WRMS | 0.00700590 m，约 7.01 mm |

SDM 正常完成并写出全部结果，但日志显示 **10000 次迭代未达到内置收敛阈值**。本表描述模板原有迭代预算、平滑和约束下的恢复；当前相关性低、零区泄漏明显，不能从良好位移拟合推断棋盘结构充分恢复。CSV 提供完整分段和数据集结果。指标基于文本输出精度：模型四位小数、预测文件单位下五位小数。
