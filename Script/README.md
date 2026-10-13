# SDM 科研脚本

这些脚本负责观测准备、棋盘试验和结果绘图，不修改 SDM 科学内核。Python 使用 NumPy；滑动及速度模型绘图另需 Matplotlib；两个 shell 绘图脚本使用 Bash、awk 和 GMT 6。

## 脚本与数据流

| 脚本 | 输入 → 输出 | 主要计算函数 |
| --- | --- | --- |
| `insar_to_sdm.py` | 方向余弦观测 → 六列 SDM 观测，无表头 | `sdm_observations` |
| `make_offset_euler.py` | 坐标或六列 SDM 观测 → 六参数修正矩阵，一行表头 | `offset_euler_matrix`、`insar_offset_euler_matrix` |
| `make_dgrn_model.py` | 用户指定的 CRUST1.0／USTC2.0 和经纬度 → DGRN2020 模型段及速度剖面图 | `crust1_to_dgrn`、`bilinear_profile`、`brocher_density`、`plot_velocity_model` |
| `make_sdm_checkerboard.py` | SDM 滑动模型 → 同几何的棋盘真值 | `make_checkerboard` |
| `run_sdm_checkerboard.py` | 反演模板和棋盘真值 → 正演、合成观测、反演、指标 CSV | `write_synthetic`、`resolution_metrics`、`fit_metrics` |
| `plot_slip.py` | SDM 滑动模型 → 3D 图、分段平面图、角点表、迹线 | `corners`、`plot_3d`、`plot_plane` |
| `plot_gps_fit.sh` | SDM 的 E/N 和／或 U 拟合输出 → 水平／垂向图、合并表 | `read_component`、`vertical_vectors`、`plot_gps_map` |
| `plot_insar_fit.sh` | SDM 的一维拟合输出 → 观测、拟合、残差三面板图 | `read_insar_fit` |

Python 脚本按读取、计算、写出组织，`main()` 负责命令行。棋盘运行器依次检查输入、正演、生成合成观测、反演、报告。`run_sdm_checkerboard.py` 复用同目录下生成器的模型读取函数，其余脚本可单独运行。shell 脚本顶部集中放置地图范围、比例、单位标签和断层覆盖层路径，数据处理函数与 GMT 绘图分开。

## 文件列与科学约定

**六列 SDM 观测与六列 SDM 拟合输出的列含义不同。** 两个 shell 脚本读取拟合输出，不能直接读取转换器生成的观测文件。

| 文件 | 列顺序 | 表头 |
| --- | --- | --- |
| 转换器原始输入 | `lon lat displacement c_up c_north c_east sigma` | 默认无；`--skiprows` 指定 |
| 固定方向 SDM 观测 | `lat lon displacement sigma` | 由模板 `nheader` 指定 |
| 逐点方向 SDM 观测 | `lat lon displacement sigma incidence azimuth` | 转换器输出无表头；模板应设 `nheader=0` |
| SDM 拟合输出 | `lat lon observation residual prediction correction` | 一行 |
| SDM 滑动模型 | `lat lon depth x_local y_local length width slip_strike slip_down_dip amplitude strike dip rake [...]` | 原生输出一行；空行分段 |

- 经纬度和角度为度；深度、局部坐标及 patch 尺寸为 km；滑动为 m。观测、误差和拟合保留输入位移单位，模板 `datunit` 表示每个输入单位对应多少米。shell 的 `UNIT` 仅改变标签。
- 入射角从上方向量起算，范围 0～180°；方位角从北顺时针起算。投影为 `d = U*cos(i) + N*sin(i)*cos(a) + E*sin(i)*sin(a)`，支持 LOS 和方位向观测；位移正负不自动反转。
- GPS 修正矩阵的行序为 **全部 N、全部 E、全部 U**，不是 ENU。参数列为 `b_N b_E b_U omega_X omega_Y omega_Z`。矩阵采用半径 6371000 m 的球体，旋转为 `omega × r`；偏移是各测点自身 NEU 中的常数。
- 修正矩阵必须与反演模板的数据集、站点和参数顺序对应。矩阵单位为 m，观测单位由 `datunit` 控制，不能直接将矩阵数值当成 cm。
- 棋盘按每段局部 patch 中心分箱，不移动或拆分 patch。SDM 输出约定为 `slip_strike = S*cos(rake)`、`slip_down_dip = -S*sin(rake)`；生成器只改滑动与 rake，应力列须经正演重算。
- 拟合为第 5＋6 列，绘图残差为第 3 列减拟合。原第 4 列不能当成观测误差。
- 棋盘运行器关闭两阶段 correction，保留模板几何、方向、误差、权重和平滑。噪声按模板数据集顺序生成；`--noise-scale 0` 不消耗随机数。详细指标与假设见 [checkerboard_guide_zh.md](checkerboard_guide_zh.md)。

`plot_slip.py` 按模型表头名称读取绘图列，棋盘生成器按 SDM 固定列位置读取并保留额外列。这两种读取方式分别服务绘图和保留原模型文本。滑动绘图的局部球面近似沿用半径 6371.0088 km；该常数与修正矩阵原有的半径分别保留。

## 使用示例

以下命令从仓库根目录运行；`raw_insar.txt`、`insar_sdm.dat` 等是示例文件名，需替换成自己的数据。

```bash
python Script/insar_to_sdm.py raw_insar.txt insar_sdm.dat
# 若原始方向余弦为 ENU 且有两行表头：
python Script/insar_to_sdm.py raw_enu.txt insar_sdm.dat -E 4 -N 5 -U 6 --skiprows 2
python Script/make_offset_euler.py insar_sdm.dat offset_insar.grn --mode insar
python Script/make_offset_euler.py SDM2025/InputFile/GPS/ns.dat offset_gps.grn --skiprows 1

python Script/make_sdm_checkerboard.py slip_model.dat checkerboard.dat \
    --cell-size 9 9 --high-slip 1 --low-slip 0 --rake 0
python Script/run_sdm_checkerboard.py template.inp checkerboard.dat --noise-scale 1 --seed 2026

python Script/plot_slip.py slip_model.dat -o figures/slip
bash Script/plot_gps_fit.sh -e east_fit.dat -n north_fit.dat -u up_fit.dat -o figures/gps
bash Script/plot_gps_fit.sh -u up_fit.dat -o figures/gps_up
bash Script/plot_insar_fit.sh -i insar_fit.dat -o figures/insar
```

棋盘运行器从**模板目录**解析模板里的相对路径，在该目录运行 SDM 并保存结果。其他脚本的相对路径从启动命令的工作目录解析。模板须指定已有的原反演滑动模型，供棋盘运行器核对几何；运行器需要已有 SDM2025 执行文件。

模板支持以 `#` 开头的注释及空行，科学记录按 UTF-8 读取；仓库原始模板在注释及被跳过的 correction 边界行中的旧编码度符号不影响科学字段的解析。

绘图前按实际区域、位移单位调整 shell 顶部配置。默认地图范围和单位沿用 Maduo 设置，仓库中的 Tohoku 样例需要相应修改。默认覆盖层为工作目录里的 `slip_model_patches.dat` 和 `slip_model_trace.gmt`；先用 `plot_slip.py` 生成，或在配置中将覆盖层路径设为空字符串。

绘图 `-o` 均为不带 `.pdf/.png` 的前缀：

| 脚本 | `-o result` 的输出 |
| --- | --- |
| Slip | `result_3d.pdf`、`result_s01.pdf` 等；`--format png` 切换格式 |
| GPS | `result_horizontal.pdf` 和／或 `result_vertical.pdf`、`result_merged.dat` |
| InSAR | `result.pdf` |

GPS 的 `-e/-n` 须成对；仅 U 也可绘图。未指定分量选项时使用配置中的默认文件；指定后只处理提供的分量。合并表的缺失分量写为 `NaN`。正 U 箭头向上、负 U 向下、零值不画箭头。Slip 所有分段共用英寸／公里比例与全模型色标，3D 横轴为经纬度、深度向下增加。

所有 Python 脚本支持 `--help`，两个 shell 脚本支持 `-h`。调试数值时可直接导入计算函数、检查输入数组和返回值；运行 SDM 出错时查看模板目录中的 `work_forward.log` 或 `invwork_inverse.log`。重跑会删除并重写本次试验的生成文件，比较多个试验时使用独立目录。

## 为 DGRN2020 提取速度模型

依赖 NumPy 和 Matplotlib。四个位置参数依次为 `mode lat lon output`，纬度在前、经度在后，单位为度，经度归一化到 `[-180, 180)`。两种模式都用 `--model` 指定输入：CRUST1 为包含 `crust1.bnds`、`crust1.vp`、`crust1.vs`、`crust1.rho` 的目录；USTC 为采用地表参考深度的五列模型文件。所有相对路径均从启动目录解析。

```bash
python3 Script/make_dgrn_model.py crust1.0 34.5 98.5 crust1_dgrn.dat \
    --model Crust/crust1.0
python3 Script/make_dgrn_model.py ustc2.0 34.25 98.25 ustc2_dgrn.dat \
    --model Crust/USTClitho2.0/USTClitho2.0.txt
```

输出含 `#` 注释、模型数据行数及五列 `no depth[km] vp[km/s] vs[km/s] rho[g/cm^3]`。将它替换到 `DGRN2020/InputFile/dgrn2020-example.inp` 的 `GLOBAL MODEL PARAMETERS` 模型段，即原行数 `7` 及后面的模型记录；保留前面的观测深度、距离／源深度采样、积分精度和输出路径，并按研究需要配置。这个 `.dat` 是模型段，不能作为完整 `.inp` 直接交给 DGRN。注释全部从首列 `#` 开始，无空行，与 `skipdoc.f` 兼容；文本使用 UTF-8 以记录用户提供的路径，数值记录仍为 ASCII。注释记录实际输入文件的绝对路径；模型或图片输出均禁止覆盖用户提供的输入模型文件。

同时在输出目录保存 `<output_stem>_velocity.png`：上述例子分别生成 `crust1_dgrn_velocity.png` 和 `ustc2_dgrn_velocity.png`。横轴为速度（km/s），红线为 Vp、蓝线为 Vs；纵轴为相对当地地表的深度（km），从顶部的 0 向下增加。图保留界面跳变，显示至底部半空间的上边界，不人为延长深度范围。梯度段按 DGRN 的连续目标剖面绘制：先线性插值 `rho`、`rho*Vp^2` 和 `rho*Vs^2`，再还原 Vp、Vs；曲线采样仅用于绘图，不参与 DGRN 的子层离散或数值输出。图片直接保存，不打开 GUI 窗口。

**正确性关键约定：**

- [CRUST1.0 官方模型](https://igppweb.ucsd.edu/~gabi/crust1.html)按 1° 单元取值，不对各层界面做空间插值。零基索引为 `ilat = floor(90-lat)`、`ilon = floor(180+lon)`，沿用 `getCN1point.f` 在整数纬度／经度边界取南侧／东侧单元的行为；南极边界归入最南单元，180° 经线归入 −180° 单元。输出注释记录实际单元中心。
- `crust1.bnds` 为海平面参考的界面海拔 `b_i`，单位 km；向下深度为 `z_i = b_0 - b_i`，层厚为 `b_i - b_(i+1)`。保留有正厚度的水以外各层（包括实际存在的冰层）、原始 Vp/Vs/密度及 Moho 以下地幔。零厚度层无需参数；速度突变界面在同一深度输出上下两侧值，连续界面只写一条。非零厚度水层会报错，避免将流体写入固体内核或擅自改变自由表面。
- [USTClitho2.0 作者说明](https://github.com/ShouchengHan/USTClitho2.0)区分地表参考和海平面参考模型。用户应通过 `--model` 指定地表参考的五列文本（`lon lat depth Vp Vs`），例如本地的 `USTClitho2.0.txt`，以适配 DGRN 的当地自由表面 `z=0`；文件名不能确定深度参考，需要用户核对。这个本地文本实际范围为经度 72–136°、纬度 18–54°，水平间隔 0.5°，深度节点为 `0, 5, 10, 15, 20, 30, 40, 60, 80, 100, 120, 150 km`；脚本读取用户文件的实际网格和深度，不补造缺失节点。每个深度按经纬度双线性插值 Vp、Vs，禁止域外外推；末节点的参数向下作为常数半空间。不会人为加入 Moho 突变。
- USTC 的密度按 [Brocher (2005), 式 (1)](https://pubs.usgs.gov/publication/70027956)计算：`rho = 1.6612*Vp - 0.4721*Vp^2 + 0.0671*Vp^3 - 0.0043*Vp^4 + 0.000106*Vp^5`，Vp 用 km/s，rho 用 g/cm³。[SCEC 的实现](https://strike.scec.org/scecpedia/UCVM_Density_Formula)给出相同公式；[USGS 的关系说明](https://pubs.usgs.gov/sir/2023/5061/sir20235061.pdf)给出 1.5–8.5 km/s 的速度范围。先插值 Vp，再计算密度。公式源于地壳岩石经验拟合，应用于上地幔节点也是显式的经验近似，不是 USTC 提供的实测密度。
- 从代码可确定：`dgmain.f` 将这些单位换算为 SI，并计算 `mu=rho*Vs^2`、`lambda=rho*(Vp^2-2*Vs^2)`。不同深度记录之间，DGRN 实际线性插值的是 **lambda、mu 和 rho**，再由 `dgsublay.f` 按 `dgalloc.f` 的 5% 参数变化尺度拆成均匀子层；并非直接线性插值 Vp/Vs。脚本检查正密度、正速度、`lambda>0`、深度顺序，并拒绝 DGRN 不接受的完全均匀模型。

数据读取和水平插值属于模型处理，重复界面记录、末行半空间及五列输出属于 DGRN 输入约定。将局部三维模型剖面用于 DGRN，是横向均匀、各向同性的一维分层弹性半空间近似；深度和断层几何需采用同一当地地表参考。界面、弹性模量、经验密度及底部半空间会影响静态差分 Green's functions，进而影响 SDM 位移预测与滑动反演。

**可复现性及验证：** CRUST1 输出保留有效层原参数和界面厚度；USTC 输出包含水平插值、密度估算及底部半空间延拓。正式反演前应核对输入文件、坐标与深度参考，并评估经验密度和底部半空间的影响。目录输入方式不改变模型数值或绘图行为。

```bash
python3 -m unittest discover -s Script -p 'test_make_dgrn_model.py' -v
```

测试覆盖界面深度／厚度、连续与突变界面、零厚度层、冰／水层、原 Fortran 单元选择、双线性插值、实际 USTC 深度、密度单位／范围、输出格式及绘图速度的模量梯度还原。
