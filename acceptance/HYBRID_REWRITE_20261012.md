# BabelMesh 混合网格重构交付与验收

复核日期：2026-10-12。分支：`codex/hybrid-remesh`。已实现新生成核心、CLI、Qt 6 CAD/网格界面、交换格式、固定域边界及基于论文的受约束质量改进。功能执行、面质量通过和 CFD 物理验证分开记录。

## 实现和调用

`schema_version: 2` 的实际顺序为：完整正交笛卡尔背景及局部加密 → 挖除实体/域外并预留贴体和层带 → 提交指定高度边界层 → 受约束空腔填充 → 唯一共享面、闭壳和体积检查 → 可选优化 → 双向采样 CAD 距离和输出。保留核心的节点不移动。

二维统一挤出一层三维体网格：六面体、三棱柱及加密接口多边形柱；front/back 为 symmetry，不能出现四面体、四棱锥或第二层厚度。三维采用笛卡尔六面体、接口四棱锥和贴体/层带四面体。支持 STEP、基本体、草图/拉伸、布尔、变换、旋转成形、放样与物理组。Qt 与 CLI 使用同一核心，任务使用工程快照，可取消。

```bash
cmake -S . -B build/qt -DCMAKE_BUILD_TYPE=Release
cmake --build build/qt -j4
build/qt/bin/babelmesh-gui examples/cylinder2d.json
build/qt/bin/babelmesh generate examples/cylinder2d.json --out outputs/cylinder2d
build/qt/bin/babelmesh check outputs/cylinder2d/mesh.json
build/qt/bin/babelmesh optimize outputs/cylinder2d/mesh.json --out outputs/optimized
build/qt/bin/babelmesh convert outputs/cylinder2d/mesh.json --to fluent-ascii --out outputs/grid.msh
build/qt/bin/babelmesh convert outputs/cylinder2d/mesh.json --to cgns --out outputs/grid.cgns
```

已有输出默认保护，需要覆盖时明确加 `--overwrite`。依赖、单位、参数和格式边界见 [README](../README.md)、[参数说明](../docs/混合网格参数与格式说明.md) 和 [重构方案](../docs/混合网格重构方案.md)。已验证 Linux/WSL + Qt 6，没有交付 Windows 原生安装包。

背景支持分段数或间距及 quadtree/octree 局部加密，满足 2:1 邻接。尺寸源分别作用于背景、表面、过渡区或全部。首层厚度、层数和增长率按壁面组设置；首层厚度不是首排单元中心壁面距离。CAD 容差独立使用米；尺寸下限、单元数、内存和线程预算独立约束。

## 固定计算域和优化

生成的 CAD 流体域必须完全位于声明背景内。二维截面加减半厚度必须落在背景及外流计算域中。检查全部输出顶点是否位于固定闭区间，不借几何容差放宽范围；每个优化候选使用同一检查，越界回滚。拒绝超出背景的 CAD 域和过厚挤出，不自动扩展背景。

Fluent/CGNS 的 BabelMesh 元数据保留锁定域；显式单位缩放同步缩放坐标、厚度和锁定范围。外部文件缺失声明边界时，独立优化锁定其初始点包围盒；不能恢复源软件未保存的原始域声明。

优化固定背景、物理边界和已标记层带节点，二维点成对移动。保留邻接、外心与面中心对齐目标及局部受约束移动。新增默认三维过渡区 3→2、4→4 内部换边；保持空腔外表面、体积、共享面方向和节点，之后按新邻接再平滑。四项全局 max/p99 面指标不得恶化，目标必须下降。实际新增节点和单元均为 **0**。可选凸多面体合并默认关闭，不能合并层带或背景。

论文原文、对应实现和限制见 [高质量混合网格论文与实现](../docs/高质量混合网格论文与实现.md)。新增区域面质量及四面体均比诊断，均比不替代 CFD 面门槛。

原生 `.bsm` 缺少区域保护标签。严格几何识别可恢复轴对齐背景和符合 front/back symmetry 的单层二维规则，但不能恢复层带高度。缺少区域标签时保留全部节点和单元，只检查和报告；继续受约束节点优化应使用完整 `mesh.json`。

## 最终验证证据

本轮最终证据在 [hybrid-20261012](evidence/hybrid-20261012)，早期和历史结果未代替本轮验收。

- CMake CLI/Qt 构建成功：[构建日志](evidence/hybrid-20261012/quality-final-build.log)。
- CTest **1/1 通过**：[契约日志](evidence/hybrid-20261012/quality-final-ctest.log)。覆盖混合体积、非法壳、二维规则、凸合并、固定边界、真实 3→2/4→4 空腔、层保护、均比归一化及 BabelSim reader。
- 集成测试 **105 次真实命令检查通过，14 个生成案例，275.70 秒**：[机器报告](evidence/hybrid-20261012/final-v2/integration_report.json)、[运行日志](evidence/hybrid-20261012/final-integration-v2.log)。包含预期拒绝的负例，不表示全部案例通过面门槛。
- 编译自实际相邻 BabelSim 的 reader 核对单元、体积及 max/p99 面几何，面指标差异不超过 `1e-8`。所用源码及可执行文件哈希分别记录。
- STEP 导出后重新导入，CAD 体积相对差异小于 `1e-10`。保留 CAD 回归还覆盖 STEP 单位、装配和 m/mm 归一化。
- Fluent ASCII、CGNS 标准混合与 NGON/NFACE 往返保持类型、连接计数、体积、物理组及自定义固定域。独立 `cgnscheck` 无 ERROR/WARNING；另用独立一四面体 Fluent 文件验证方向与体积 `1/6`。显式缩放核对体积三次方关系与边界同步变化。
- 负例覆盖越界 CAD、过厚挤出、越界检查点和优化、二维四面体、第二厚度层、层间隙、资源预算、严格面质量、Fluent HDF5、覆盖保护与路径别名。
- Qt 实际编辑、保存、后台生成得到 577 个单元：[保存工程](evidence/hybrid-20261012/gui-release/workflow/edited-project.json)、[网格报告](evidence/hybrid-20261012/gui-release/workflow/generated/mesh_report.json)、[截图](evidence/hybrid-20261012/gui-release/workflow.png)。工作流改了表面尺寸，与原始 580 单元案例不同。
- Qt 真实任务取消：[报告](evidence/hybrid-20261012/gui-release/cancel/cancellation-report.json)，进程停止且未发布结果。
- 五个最终 GUI 捕获覆盖 CAD、二维工作流/等轴透视、三维正交/透视：[汇总](evidence/hybrid-20261012/gui-release/gui_report.json)。检查实际 framebuffer 的颜色、透明合成、黑色边及剖切共面；采用 Xvfb/软件 OpenGL，没有自动逐一点击所有 CAD 按钮。
- 旧版 `make MODE=release test -j4` 回归通过：[日志](evidence/hybrid-20261012/legacy-final-regression.log)。最新 CLI 的 Make 构建与 M0 合同另有 [日志](evidence/hybrid-20261012/make-final-cli.log)。
- 实际暂存源码导出到独立目录后重新 CMake 编译，[构建](evidence/hybrid-20261012/delivery-build.log)、[CTest](evidence/hybrid-20261012/delivery-ctest.log)、[CLI/reader/CGNS 运行](evidence/hybrid-20261012/delivery-runtime.log) 及 [独立源码 M0 13 项合同](evidence/hybrid-20261012/delivery-m0.log) 均通过。仍需安装文档指定的外部开发库和测试用 BabelSim 源码；不声称二进制或系统库自包含。

## 最终网格质量

面门槛不变：max/p99 非正交角 **75°/60°**，max/p99 内部归一化偏斜 **0.5/0.25**。以下通过闭壳、正体积、共享面、声明边界、固定域、双向采样 CAD 距离和实际 reader。面门槛单独列出。

| 案例 | 单元数 | 单元组成 | 非正交 max/p99 (°) | 偏斜 max/p99 | 面门槛 |
|---|---:|---|---|---|---|
| 二维圆柱，两层 | 580 | 340 六面体 + 240 三棱柱 | 57.484 / 57.484 | 0.2068 / 0.2015 | 通过 |
| 二维局部加密 | 688 | 436 六面体 + 12 多面体 + 240 三棱柱 | 57.484 / 57.484 | 0.2455 / 0.2125 | 通过 |
| 三维内流盒子 | 216 | 216 六面体 | 0.000 / 0.000 | 0.0000 / 0.0000 | 通过 |
| 三维球体 | 2,381 | 920 六面体 + 128 四棱锥 + 1333 四面体 | 68.518 / 49.881 | 0.4258 / 0.2664 | p99 偏斜未通过 |
| 三维球体，两层 | 4,711 | 841 六面体 + 192 四棱锥 + 3678 四面体 | 74.631 / 71.209 | 0.4698 / 0.2867 | p99 非正交与偏斜未通过 |
| 三维球体，局部加密 | 15,978 | 720 六面体 + 774 四棱锥 + 14484 四面体 | 66.602 / 49.941 | 0.4602 / 0.2954 | p99 偏斜未通过 |

x/y 截面二维内流通道各 64 个六面体、布尔孔洞二维案例 96 个单元也通过本轮面门槛。三维锐盖圆柱、穿盒圆柱和 STEP 球体仍是诊断结果。完整质量、层高及 CAD 误差见 [集成报告](evidence/hybrid-20261012/final-v2/integration_report.json)。

二维两层实测法向厚度 0.04/0.048 m，总厚 0.088 m，增长率 1.2，覆盖 1，误差接近浮点舍入。三维首层目标 0.035 m，实测约 0.034850–0.035156 m，p95 相对误差约 0.2553%，两层总厚 0.077 m，覆盖 1。层高通过与面质量通过是两个不同结论。

与上一轮相同生成参数（仅输出目录不同）的比较见 [机器对比](evidence/hybrid-20261012/quality-comparison.json)。旧基线来自完整通过的 99 命令记录；原始质量数据另保留在 [基线报告](evidence/hybrid-20261012/verified/integration_report.json)。

| 案例 | 单元数 | 最大偏斜 | p99 偏斜 | p99 非正交角 |
|---|---:|---|---|---|
| 三维球体 | 2,396 → 2,381 | 0.4356 → 0.4258 | 0.2777 → 0.2664 | 51.345° → 49.881° |
| 三维球体，两层 | 4,722 → 4,711 | 0.4698 → 0.4698 | 0.2996 → 0.2867 | 71.209° → 71.209° |
| 三维球体，局部加密 | 16,023 → 15,978 | 0.5507 → 0.4602 | 0.3034 → 0.2954 | 50.343° → 49.941° |

薄层案例的高 p99 非正交角位于固定四面体层带。过渡区改善后它仍未下降；不能用减少单元或其他指标改善来替代该项通过。严格模式继续拒绝相关结果。

## GUI 显示

磨砂半透明材质、柔和光照、不透明度、正交/透视和黑色实线边已启用；透明合成修正了原来偏色的问题。三维使用真实平面剖切及共面截面，支持 X/Y/Z 轴和整单元筛选，取代原来按单元中心截断造成的锯齿形截口。剖切时不收缩显示单元；完整视图的显示间隙不改变源坐标。隐藏面边通过深度检查遮挡，三维实体壁面用中性色与过渡单元区分。

[三维透视截图](evidence/hybrid-20261012/gui-release/sphere3d-perspective.png) 和 [三维正交截图](evidence/hybrid-20261012/gui-release/sphere3d-orthographic.png) 来自最终实际网格。

## 格式与尚未完成的范围

BabelSim v3 面环及 owner/neighbour/patch 足以承载这些混合单元，不需增加独立二维文件类型。支持 checkpoint/BabelSim v3/Fluent ASCII 双向、CGNS 单 base 单非结构 3D zone 的标准一阶混合或 NGON/NFACE 双向、VTU 输出。只转换体网格、坐标和边界组，不转换场、材料、求解模型或数值设置。

没有实现 `.ccm/.sim/.hex/.igg` 或 Fluent binary/HDF5 原生适配，也没有启动商业软件的最新版本做导入验收。STAR-CCM+、FINE/Marine 须导出受支持的交换体网格；不能把某个产品的曲面 CGNS 接口当成任意体网格接口。

尚未覆盖复杂层终止/相交和尖角截断、全局单元交叠扫描、多区/结构/overset/周期 CGNS、高阶单元、完整 CFD profile 和实际流动精度验证。`full_cfd_profile_evaluated: false` 保留。CAD 距离是双向采样到精确 B-Rep 的检查，不是连续误差上界证明。历史 `verify-suite` 仍为预留命令，新线使用 `make test-hybrid`。

交付包含源代码、文档、可复现输入和紧凑最终证据，不包含批量大网格或构建缓存。哈希见 [源码清单](evidence/hybrid-20261012/source-manifest.json) 和 [二进制哈希](evidence/hybrid-20261012/binaries.sha256)。原始工作区基线归档保留在 `/home/midway/.cache/babelmesh-baselines/20261009/before-hybrid-rewrite.tar.gz`，未纳入提交。
