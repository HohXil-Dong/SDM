# SDM 科研脚本

这些脚本负责观测准备、棋盘试验和结果绘图，不修改 SDM 科学内核。Python 使用 NumPy；滑动绘图另需 Matplotlib；两个 shell 绘图脚本使用 Bash、awk 和 GMT 6。

## 脚本与数据流

| 脚本 | 输入 → 输出 | 主要计算函数 |
| --- | --- | --- |
| `insar_to_sdm.py` | 方向余弦观测 → 六列 SDM 观测，无表头 | `sdm_observations` |
| `make_offset_euler.py` | 坐标或六列 SDM 观测 → 六参数修正矩阵，一行表头 | `offset_euler_matrix`、`insar_offset_euler_matrix` |
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
