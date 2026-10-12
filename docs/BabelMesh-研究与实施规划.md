# BabelMesh：三维非结构网格生成器研究与实施规划

日期：2026-09-23。本文记录研究结论、目标架构和验收门禁；实现状态以
`acceptance/suite.json` 与实际执行日志为准。已实现的功能只按对应门禁证据标记，
后续阶段的规划文字不代表功能已经存在。

## 1. 建议与边界

采用 C++17、Open CASCADE Technology（OCCT）几何内核，开发独立、可批处理的三维网格生成器。主路线为：均匀笛卡尔六面体背景网格 → 几何与区域自适应 → 流体域提取 → 特征约束贴体 → 质量优化 → 核心网格变形与边界层插入 → 质量门禁 → BabelSim v3 导出与求解验收。

目标产物是六面体占主体、允许显式多面体过渡和局部棱柱层的贴体有限体积网格。初始网格必须是真正的三维均匀 Cartesian hex；不能以四面体网格加格式转换代替该路线。粗细过渡单元即便几何外形仍像立方体，拆面后也可能超过 6 个面，必须按多面体统计，不能冒充严格全六面体。

V1 必须覆盖：基本点线面体 CAD、布尔运算、STEP 导入、内外流体域、局部区域/曲率/间隙加密、壁面边界层、方向性区域加密、自动质量报告、BabelSim 读取与数值验证。基础 CAD 首先以 C++ API 和声明式工程文件交付；GUI 是后续可选层。用户未要求完整约束草图器、装配约束求解器或专业 CAD 编辑器。

不把任意损坏 STEP 自动修好、任意尖角实现指定层数、复杂几何必然全六面体作为 V1 保证。有效输入仍须属于验收覆盖的几何范围；对超出范围或无法满足尺寸/质量约束的输入，必须自动明确失败。不得静默删除小特征、封堵真实通道、放大首层或取消壁面层后报告成功。

## 2. BabelSim 当前源码契约

本次只读核对的 HEAD：`529daf898c868ebbac64fb4ec9daafb15dbaa094`，工作区无未提交修改。没有运行编译、回归或新的数值算例。历史笔记中的 v2/仅六面体结论已被当前源码替代。

### 2.1 文件组织

`case.bs` 通过 `mesh mesh/xxx.mesh` 指向网格；物性、初场、离散格式、求解参数与输出配置另存。网格文件只描述几何、拓扑与 patch，不保存速度、压力、材料或边界条件数值。

证据：[读取器](/home/midway/BabelSim/src/io/mesh_reader.cpp:44)、[接口声明](/home/midway/BabelSim/include/babelsim/mesh_io.h:9)、[示例 case](/home/midway/BabelSim/cases/heat/case.bs:1)。

当前唯一运行时格式：

```text
BABELSIM_MESH 3
vertices N
x0 y0 z0
...
faces F
face K owner neighbour patch v0 v1 ... v(K-1)
...
patches P
patch patch_name patch_kind
...
end
```

- 所有 ID 从 0 开始，顶点、面、patch ID 由记录顺序确定；cell ID 由 owner/neighbour 引用确定。
- 没有 cells 块，也不写 cell 类型、体积和中心。单元数为最大被引用 cell ID 加 1；生成器必须保证 ID 紧凑、无孤立单元。
- 内部面：`owner >= 0`，`neighbour >= 0`，两者不同，`patch = -1`。
- 边界面：`owner >= 0`，`neighbour = -1`，`patch` 为有效 patch ID。
- 每个面只写一次，一个 owner 和至多一个 neighbour。多个子面连接同一对单元是允许的。
- K 至少为 3，面顶点环不能重复顶点。K 可变，单元面数也可变。
- patch kind 可为 `generic/wall/inlet/outlet/symmetry/processor`；reader 也把 `mirror` 读作 symmetry。生成器导出全局物理网格不创建 processor patch。
- parser 按空白分词，不支持注释，也不接受 end 后额外 token。patch 名不能含空格；生成器额外保证唯一、稳定、可追踪。
- N/F/P 必须为正。当前 Index 为 int32；导出前不仅检查顶点/面/cell 数，还要检查 CSR 总连接数、offset 和加一溢出，超限明确失败。
- 文件没有单位标记。工程文件内部统一 SI 米，导出前完成缩放，在旁文件记录源单位与变换。
- 建议浮点输出 `max_digits10`，固定 locale，稳定排序；不能为使文件更短而损失薄层坐标精度。

一格单位立方体的完整可读示意（不是本轮生成的网格产品）：

```text
BABELSIM_MESH 3
vertices 8
0 0 0
1 0 0
1 1 0
0 1 0
0 0 1
1 0 1
1 1 1
0 1 1
faces 6
face 4 0 -1 0 0 3 2 1
face 4 0 -1 0 4 5 6 7
face 4 0 -1 0 0 1 5 4
face 4 0 -1 0 1 2 6 5
face 4 0 -1 0 2 3 7 6
face 4 0 -1 0 3 0 4 7
patches 1
patch walls wall
end
```

### 2.2 内存拓扑和方向

[mesh.h](/home/midway/BabelSim/include/babelsim/mesh.h:153) 中用 CSR 存储 face→vertices、cell→faces/signs、cell→vertices、cell→neighbours；另有 owner/neighbour/patch 数组。面面积向量朝 owner 外侧，owner 的 cell-face sign 为 +1，neighbour 为 -1。几何构建计算单元体积/形心、面中心/面积向量、插值权重、非正交修正等。

读取器会调用 `Mesh::polyhedral()`，构造时自动调整面环方向。生成器仍应输出一致的正确方向，不能把自动翻面当成修复任意坏单元的手段。

[mesh.cpp](/home/midway/BabelSim/src/core/mesh.cpp:319) 的关键约束：

1. 有限坐标，非退化面，正的有符号单元体积；体积/形心来自面三角化后的有符号四面体积分。
2. 内部面必须满足 `Sf · (Cn-Co) > 0`；边界面必须满足 `Sf · (Cf-Co) > 0`，且相应向量长度非零。
3. 每个单元的有向面面积向量闭合；现源码阈值为 `1e-9 * max(sum(faceArea), 1.0)`。
4. 每个单元的闭合外壳中，每条无向边恰好出现两次且方向相反。
5. 凹多边形面有 ear-clipping 处理，但“可变多面体拓扑”不等于任意凹形、自交或病态几何均可用于 FVM。其翻面启发式和中心方向条件对严重凹单元尤其敏感。

因此本项目应优先生成凸或温和非凸单元；强凹单元拆分或回退细化。生成器的几何自交、体重叠、重复面检测应强于求解器当前 validate，不能把 reader 通过等同于几何完全有效。

### 2.3 粗细网格交界的决定性要求

当一个粗网格面邻接四个细单元时，粗面须替换为四个共享子面；细侧与粗侧引用同一组 face 对象。粗单元其余相邻面也必须补入分裂边上的中点，保证所有单元内的边段配对。

仅将一个大面拆成四面，而保留邻面的整条长边，仍会触发非流形边错误。这是 M2/M3 的独立硬验收项。

还有一个当前实现特有的限制：[凸面三角化](/home/midway/BabelSim/src/core/mesh.cpp:120) 固定从第一个顶点作扇形三角化，而 [polygonGeometry](/home/midway/BabelSim/src/core/mesh.cpp:174) 拒绝零面积三角形。因此，在面环补入共线中点后，某些起点顺序会产生退化扇形，即使多边形本身合法也会被拒绝。导出适配器必须模拟实际三角化进行预检；必要时把该平面多边形一致拆成非退化三角面。对凸平面面，可使用内部面心连向每条边段的扇形；面心为新共享顶点，内部面只拆一次。不能简单删共线点而破坏相邻单元的边段一致性。

`polygonTriangles()` 的面积判定使用 `1e-13 * max(face_extent, 1.0)^2`，含绝对尺度下限。极小物理尺度、极薄层的面可能在正面积情况下被拒绝。M0应测定目标单位/尺寸范围的实际可读性，按实际导出坐标检查；不能只用理想无量纲质量标准，或擅自缩放几何而不一致调整物理量。超出当前求解器尺度范围应报告独立兼容问题。

证据：[闭合和逐边验证](/home/midway/BabelSim/src/core/mesh.cpp:709)、[同一对单元多个子面耦合测试](/home/midway/BabelSim/tests/assembly_solver_test.cpp:12)。这些是源码/测试设计证据，本轮没有重跑该测试。

### 2.4 并行边界

[parallel_context.cpp](/home/midway/BabelSim/src/parallel/parallel_context.cpp:512) 当前由 rank 0 读全局网格，随后广播全局网格并按邻接关系构造各 rank 局部网格。BabelMesh V1 只需导出一个全局 v3 文件，不生成 MPI ghost 层。BabelMesh 的多线程扩展性与 BabelSim 的并行读取内存代价必须分别测量。

## 3. 公开算法依据与本项目选择

NUMECA 官方案例明确描述了 HEXPRESS 的 volume-to-surface、非共形全六面体、特征贴体、正交性优化以及缓冲层插入思想；Cadence 目前公开资料同时介绍带悬挂节点的全六面体、共形混合网格、inflation/extrusion 边界层及方向性加密。这里借鉴公开思想，不声称知道其闭源细节，不混同 HEXPRESS、HEXPRESS/Hybrid 和 Pointwise 的不同实现。[NUMECA 官方案例，第 2–3 页](https://www.cadence.com/content/dam/cadence-www/global/en_US/documents/tools/system-analysis/paper-generic-gas-turbine-convertor.pdf)，[Cadence 网格产品说明](https://www.cadence.com/en_US/home/tools/system-analysis/computational-fluid-dynamics/pre-processing-meshing.html)。

初始均匀 hex → 细化 → 贴体 → 优化 → 黏性层是本项目确定采用的阶段设计。2:1 平衡、具体目标函数、回退次序及本文数值阈值均为 BabelMesh 工程选择，不应标成 HEXPRESS 官方参数。

OpenFOAM 官方对 castellation、feature snapping 和 layer addition 的公开说明，可作为独立实现的补充算法参照；尤其是 seed 保留流体连通区和通过核心网格收缩腾出插层空间。引用思路不等于直接复制其代码或把调用 snappyHexMesh 作为 BabelMesh 自研完成。[Castellation](https://doc.openfoam.com/2606/tools/pre-processing/mesh/generation/snappyhexmesh/castellation/)，[Snapping](https://doc.openfoam.com/2306/tools/pre-processing/mesh/generation/snappyhexmesh/snapping/)，[Layer addition](https://doc.openfoam.com/2212/tools/pre-processing/mesh/generation/snappyhexmesh/layers/)。

严格全六面体的复杂几何拓扑处理应视为独立研究目标；六面体网格综述与 HybridOctree_Hex 可作为后续 transition template/Jacobian 优化参考。[Hex-Mesh Generation and Processing: a Survey](https://arxiv.org/abs/2202.12670)，[HybridOctree_Hex](https://arxiv.org/abs/2401.05984)。

本次检索定位到 Cadence Marine 教程，但全文重取遇到 timeout/403；不以未取得的教程细节支撑算法结论，也不依赖第三方转载手册中的性能承诺。

## 4. CAD 与几何模型

### 4.1 工具链

采用 C++17 与 GNU Make 管理编译目标、依赖检查、验收测试、Debug/Release 构建和清理。CAD 使用 OCCT；工程配置使用带 schema_version 的 JSON，数值计算使用 Eigen。版本和许可证记录在 `THIRD_PARTY_NOTICES.md`，先基于 Ubuntu 22.04 的可构建依赖固定基线。几何谓词和拓扑代码禁止 fast-math，以免破坏退化判定；与 BabelSim 可通过进程/文件集成，避免继承其编译选项。

OCCT 提供 B-Rep 建模与 STEP 能力。拟使用 BRepBuilderAPI、BRepPrimAPI、BRepAlgoAPI、BRepCheck、ShapeFix、BRepMesh、STEPCAFControl/XCAF 等模块，具体 API 以锁定版本文档和编译结果为准。[OCCT 建模算法](https://occt3d.com/dev/doc/overview/html/occt_user_guides__modeling_algos.html)，[OCCT STEP 文档](https://github.com/Open-Cascade-SAS/OCCT/blob/master/dox/user_guides/step/step.md)。

三角面 BVH、相交和最近点查询可先用 OCCT/自有封装；如需要更强谓词或面自交检测，再引入 CGAL 的具体包。精确谓词并不自动使所有浮点交点构造精确；相交构造仍须统一事件 ID、容差和必要的精确构造。引入前核对具体依赖包的许可证和发布兼容性，不直接假定与仓库许可一致。[CGAL Polygon Mesh Processing](https://doc.cgal.org/latest/Polygon_mesh_processing/index.html)。

### 4.2 必备建模操作

| 层级 | V1 能力 | 主要验收 |
|---|---|---|
| 点/线 | 3D 点、线段、圆弧/圆、插值样条、闭合 wire | 端点、闭合性、长度可核对 |
| 面 | 平面/带孔平面、由边界构面、放样或拉伸得到曲面 | 面积、法向、边界环合法 |
| 体 | box、sphere、cylinder、cone、extrude、revolve、loft；合法 shell 成 solid | 体积、闭合、连通分量可核对 |
| 变换 | 平移、旋转、统一缩放、镜像、复制 | 单位、法向、体积变化正确 |
| 布尔 | fuse、cut、common；split 用于划 patch 或分区域 | 成功状态、体积、拓扑历史可核对 |
| STEP | .step/.stp 导入、装配实例变换、名称/颜色尽力保留、单位转换、实体选择 | 不仅导入三角显示，必须得到可查询 B-Rep |
| CFD 域 | 外流外包域减固体；内流选流体体积并封闭入口/出口 | 流体 seed、边界 patch 和通道连通正确 |

导入流程：读入 → 单位/实例变换 → 合法性检查 → 在显式容差内缝合/修复 → 再检查 → 建立实体命名 → 离散用于加速查询。不得“读成功”就直接网格化。

边界标签采用操作历史＋几何选择器，不使用不稳定的 face 数组下标。布尔后的 Modified/Generated/Deleted 关系应能传播标签；出现一面多标签冲突或一对多歧义时，输出可定位诊断并要求工程文件规则解决。导出 mesh face 保存 CAD 关联到旁文件，便于追踪。

几何离散必须在共享 CAD 边上统一采样，避免每个面独立三角化留下裂缝。保留 face/edge/vertex 归属和参数坐标；光滑曲面接缝不是自动锐边，需区分周期 seam、真实棱边、patch 接缝和自由边。

## 5. 核心网格算法

### 5.1 背景网格与尺寸场

以用户计算域为边界生成均匀三维 Cartesian 网格，每轴步长恒定。默认近等边；严格立方体模式通过调整外包域实现，不能为了严格立方体偷偷改变用户固定边界。

初始 `(nx,ny,nz)` 或 `h0` 必须显式记录。基础索引使用整数格坐标；自适应采用 root-grid + sparse octree forest，稳定键由 root ID、层级、整数位置组成。树深与 Morton 位宽必须检查。

标量尺寸场建议：

`h_target(x) = min(h_global, h_surface, h_curvature, h_gap, h_region, h_feature)`。

规则覆盖 box/sphere/cylinder、距点线面距离带、命名 CAD 面组、尾流走廊及水线带。小尺寸优先，另有显式禁加密/硬上限冲突检查。曲率规则可按弦高容差推导 `h <= sqrt(8*epsilon_chord/|kappa|)`，仅作小弦长近似；最终以实际几何距离核验。间隙按相对不同表面的距离估计，不能误把同一光滑曲面附近三角形当成窄缝。

规定 N_gap（建议基础夹缝基准至少 6 个跨缝单元），根据 gap/N_gap 细化。碰到未解析薄壁/狭缝，先细化，不得仅凭 cell 中心 inside/outside 删除薄几何。区域判定使用包围盒/交叠查询，不能只采样中心而漏掉薄加密带。

M3 采用各向同性 1→8 细分，叶单元邻接层级差最多 1；边/角邻接也作保守平衡，降低复杂边环风险。尺寸场平滑与过渡宽度要独立控制；2:1 拓扑平衡不等于最终非正交性合格。

项目尺寸源按 `kind` 使用 `box`、`sphere`、`cylinder`、`cad_group`、`distance`、`curvature`、`gap`。CAD 距离与间隙通过 OCCT 的单元包围体—命名组精确距离计算，包围盒只作候选筛选；`gap` 要求两个不同 CAD 组，并要求单元同时落在两组的 `region.distance` 距离带内。`curvature` 只保留选定 CAD 面上主曲率半径不大于该源 `target_size` 的面，再按距离带加密。曲率当前按修剪面内的参数采样估算，因此后续高曲率尖点的鲁棒性仍需加密采样和曲面基准验证。所有 CAD 源只改变背景网格尺寸场，不代表已完成贴体或流体域裁切；这属于 M4。

### 5.2 拓扑组装

整数格节点去重 → 叶网格邻接查询 → 粗细界面覆盖拆分 → 相邻面边环同步细分 → 唯一面表 → owner/neighbour → cell-face CSR → 方向和闭合检查。

同一几何面不能粗侧一份、细侧四份同时保留；不能输出隐藏 hanging-node 插值关系，因为 BabelSim v3 无此字段。贴体后的点去重依赖原有拓扑/几何事件来源，不以全局宽容差将靠近的两侧薄壁焊死。

### 5.3 流体域提取

对闭合 B-Rep 和加速三角面进行单元/边界相交判断，区分 Inside/Outside/Intersect/Uncertain。使用显式 fluid seed 和被几何阻断的邻接 flood fill 确定目标连通域；入口/出口必须有盖面或等价封闭的流体体定义。

交切单元先充分细化，再形成可贴体的阶梯边界。禁止把所有 intersect 单元一删了之当作完成贴体。不能用体素分辨率吞掉薄实体，也不能自动保留所有封闭空腔。闭合性、预期通道、种子数量、域体积都进入验收。

V1 主路线为 castellated + snap，而不是 exact cut-cell。后续若引入 cut-cell，必须单独设计微小切割体合并、凹形处理与边界层连接，不能作为无验收的自动后门。

当前代码的 M4 入口为 `babelmesh generate PROJECT.json --cad-staircase`，CAD 表面缓冲带的最大单元边长由显式 `sizing.surface_size` 控制；贴体项目未配置时立即报错，不再隐式取 `minimum_size`。CAD拟合容差通过独立的 `quality.max_surface_deviation_m` 设置。实现包含 B-Rep/BVH 单元分类、尺寸加密、seed flood fill、笛卡尔阶梯抽取、最近修剪 CAD 面投影、双向采样距离审计，以及保持正体积的回退搜索。采样误差仍非连续Hausdorff上界，也不覆盖自交/重叠检查；旋转/曲面基准及完整M4必需算例通过前，M4状态维持 `in_progress`。

### 5.4 特征贴体与质量优化

建立边界点约束：CAD 顶点锁定、特征边上一维移动、CAD 面上二维投影、内部点自由。以 patch 归属和几何拓扑限定投影候选，防止薄壁两侧吸到同一表面。

分步目标位移 → 弹性/加权 Laplace 传播 → 带步长搜索的质量门禁 → 接受或回滚。先消除负体积/自交和方向错误，再改善非正交、skew、面翘曲与尺寸变化；优化不能靠牺牲几何偏差或第一层高度达标。

锐角、尾缘、patch 交线需要捕捉节点与局部面拓扑调整，不是只移动已有点就能普遍解决。阶段验收要包含旋转方盒、圆柱、圆角、薄间隙、有限厚度尖尾缘。未成功的点/面 ID、CAD ID、误差及重试原因进入报告。

高度翘曲多边形可以一致拆成平面三角面；内部共享面必须全局只拆一次，让两个单元复用相同三角化。拆面不会自动解决单元重叠或负体积，随后必须重新检查。

### 5.5 边界层：必须处理核心空间

边界层不是“在已有贴体网格上直接多放一些单元”。先把原壁面到总厚度 T 的空间从核心网格中腾出来：壁面几何仍固定，核心边界向流体内部移至层顶面，内部顶点随之变形，然后在真实壁面与层顶面之间构造层柱。外流壁面朝物体外部生长，内流壁面朝管腔内部生长；不能直接使用 solid 的外法向而不检查流体侧。

用户指定首层完整厚度 `t1`、层数 N、增长率 r，或者以明确的组合求剩余量：

`t_k = t1 * r^(k-1)`；`T = t1*(r^N-1)/(r-1)`，r=1 时 T=N*t1。

首个 cell 中心壁距与完整层厚不同，正交平面上约为 t1/2。如果配置用 y+，由 `y_center = y+ * nu/u_tau` 估计，再推完整首层厚；必须记录 u_tau 的输入或估计依据，并在流场计算后核查真实 y+。网格器不能凭目标 y+ 保证最终 y+。

实现次序：

1. M5：平壁、平行通道和光滑圆柱壁，统一层数，无局部终止。壁面及 patch 接缝协调位移，层柱各层共享拓扑，层顶与核心同一面表。四边形基面得到 hex 层，三角形基面得到 prism 层。
2. M6：曲面法向平滑、曲率约束、对置壁碰撞、局部层厚上限。通过 BVH 面相交和位移路径检查，不只测最终节点距离。
3. 在允许的策略内，优先保持 t1，减少层数/总厚度；局部减少必须通过终止楔形/多面体过渡模板并同步相邻层柱，禁止突然删除层柱留下裂口。
4. strict 模式指定区域必须达到层数/首层/覆盖率，否则失败；permissive 模式可输出降级候选与详单，但不能作为 strict 验收通过。
5. 窄缝预留核心通道，例如 `T_left + T_right + h_core_min <= gap`；这只是预检查，最终仍核验真实几何距离、层柱相交和最小通流截面。

边界层外过渡必须比较层顶法向厚度与核心法向尺度，不能拿长边尺寸来冒充匹配；过大的跳变触发核心加密/过渡重建。

### 5.6 海洋应用的方向性区域加密

严格 octree 的 1→8 只能等比例细化。M7 扩展到带每轴层级的 Cartesian refinement tree，支持 x/y/z 轴方向二分；每个轴相邻尺寸比受限，face-plane overlay 生成通用矩形子面并同步边环。

用于水线薄带、尾流、船体纵向/法向不同尺度。早期先实现轴对齐方向性；任意旋转度量张量网格不纳入 V1。方向性细化后已不是单一层级 octree，接口与统计不能继续假称各向同性。

## 6. 软件分层与产物

```text
CLI / C++ API / JSON project
  -> CAD model + repair + physical groups
  -> geometry queries + feature graph
  -> sizing field + Cartesian refinement tree
  -> shared face topology + fluid region selection
  -> feature snapping + core quality optimization
  -> layer planning + core displacement + layer topology
  -> independent quality validator
  -> BabelSim v3 + VTK + provenance / reports
```

建议模块：`geometry`、`sizing`、`cartesian`、`topology`、`snapping`、`layers`、`quality`、`io`、`cli`、`integration`。核心库不依赖 Qt、BabelSim 私有数组或 MPI 求解运行时。集成测试可链接 BabelSim 公共 Mesh API，正式生产流程通过文件交接。

内部使用 64 位 ID，最终 compact/freeze 成连续数组；v3 writer 检查所有 32 位边界。CAD 关联、层柱、细化树等开发数据不塞进 v3。

拟交付命令（目前未实现）：

```bash
babelmesh inspect-cad examples/step_external/project.json --out outputs/inspect
babelmesh generate examples/step_external/project.json --out outputs/step_external
babelmesh check outputs/step_external/mesh.bmesh --profile strict
babelmesh export outputs/step_external/mesh.bmesh --format babelsim-v3 --out outputs/step_external/mesh.mesh
babelmesh verify-suite acceptance/suite.json --out outputs/acceptance
```

`mesh.bmesh` 是规划中的带版本内部 checkpoint，不是 BabelSim 读取格式。M0 固化 schema 和 CLI，后续不得漂移。

输出至少包含：`mesh.mesh`、`mesh.vtu`、各阶段 checkpoint、`geometry_report.json`、`mesh_quality.json`、`boundary_layers.json`、`boundary_map.json`、`manifest.json`、日志和失败实体可视化。manifest 包含配置/输入 hash、依赖版本、单位、提交、线程数、时间/内存、重试和最终状态。vtu 支持真实多面体，不将所有 cells 标成 hex。

失败返回非零，不覆盖已有通过验收的结果；可把失败候选放在独立 failed 目录。相同配置和固定种子单线程要求稳定拓扑/文件 hash，多线程先要求规范化拓扑一致与几何误差有界，不在无证据时承诺跨平台逐位相同。

## 7. 质量定义与建议验收阈值

以下为本项目拟定的起始门禁，不是 HEXPRESS 指标或所有 CFD 的通用合格线。应在 M0 形成机器可读 profile；后续变更必须说明原因、影响基准和前后证据，禁止为通过验收临时放宽。

设 Sf 为有向面积向量，Co/Cn 为相邻单元体积形心，Cf 为面中心。内部 d=Cn-Co；边界 d=Cf-Co。

- 非正交角：`theta=acos(clamp(Sf·d/(|Sf||d|),-1,1))`。
- 内部面 skew：直线 Co+t*d 与过 Cf、法向 Sf 的平面相交于 P；`s=|Cf-P|/|d|`。边界 skew 另报投影位移/局部切向尺度，不能沿用边界 d 导致恒为零。
- 闭合：`e_close=|sum(sign*Sf)|/sum(|Sf|)`；生成器使用无量纲阈值，并另执行 BabelSim 原验证。
- 体积比：`max(Vo,Vn)/min(Vo,Vn)`，按背景、过渡、层区分别报告。
- 几何误差：双向表面距离的自适应采样最大值及采样分辨率；不能把仅检查网格节点称为严格 Hausdorff 上界。

| 指标 | 起始门禁 |
|---|---|
| 负体积、零面积、NaN/Inf、重复面、重叠体、非流形/开裂 | 0 个；不能用平均值掩盖 |
| 所有 owner/neighbour 投影 | 严格 >0，并报告归一化最小值 |
| 单元面积闭合 e_close | <=1e-10；另通过原生 reader/validate |
| 非正交 | max<=75°，按面数 p99<=60°；内部/壁面/过渡各报 |
| 上述定义的内部面 skew | max<=0.5，p99<=0.25 |
| 相邻体积比 | 普通各向同性核心/过渡建议<=12；薄层按法向增长独立门禁 |
| 规则 hex 的 scaled Jacobian | >0，基准建议 min>=0.1；不能套给任意多面体 |
| 曲面几何偏差 | 基准 eps_geom<=min(1e-4*L,0.05*h_surface)；另细查窄缝与局部特征 |
| 基础 BL 首层厚度 | 面积加权 p95 相对误差<=10%，max<=20% |
| 基础 BL 增长 | 指定 r 常用1.1–1.25，实测最大<=1.3；r=1允许 |
| 基础 BL 完整 N 层覆盖 | 指定光滑壁面面积100% |
| 复杂 BL 覆盖 | 可生成至少一层面积>=99%，完整请求 N 层面积>=95%；strict指定区仍100% |
| patch 覆盖 | 每个边界面恰有一个合法标签；目标 patch 不得漏失/歧义 |

高长宽比是边界层的预期性质，不统一以 aspect ratio<10 限制。对最小体积采用局部尺度与机器精度相适应的阈值，防止把不同单位或合法薄层误判；同时报告原始体积和局部归一化体积。

必须分别报告按面数的质量分位数、按面积的边界覆盖、按单元数和体积的严格 hex 比例，保留最差实体 ID/空间位置。不得用“99%单元很好”掩盖一个负体积。

## 8. 验收矩阵与阶段门禁

| 阶段 | 开发交付 | 必须通过的证据 |
|---|---|---|
| M0 契约 | 源码接口核对、格式/schema/CLI、质量profile、基准定义 | 现 HEAD 与源码定位；基准输入/期望/误差定义完整 |
| M1 CAD | 点线面体、变换、布尔、STEP、域选择、标签 | 解析体积/面积；STEP单位与装配实例；无效壳拒绝；布尔标签追踪 |
| M2 基础网格 | 均匀Cartesian、CSR、v3/VTK、独立checker | 1格/2格/多格解析体积；内部面只写一次；真实BabelSim读取；故意错误文件拒绝 |
| M3 局部加密 | octree、2:1、尺寸源、粗细共享子面 | 1粗对4细、边/角细化、薄带穿过cell角点；体积守恒；无裂缝；reader通过 |
| M4 贴体 | 几何关联、流体提取、特征捕捉、优化/回滚 | 旋转box、sphere、cylinder、圆角、窄缝、有限厚尾缘；几何偏差与质量门禁 |
| M5 基础BL | 总厚度规划、核心变形、统一层柱 | 平行板/圆柱，10层、r=1.2、完整覆盖；厚度/法向/体积/核心接口正确 |
| M6 复杂BL | 凹凸角、碰撞、逐层终止模板、局部回退 | 对置窄缝、凹角、三面交角；无穿层；覆盖率；不满足条件明确失败 |
| M7 方向性加密 | 每轴细分、各向异性界面、自由面/尾流带 | z加密不无故等量细分x/y；接口闭合；与各向同性结果对比单元数和质量 |
| M8 集成发布 | 全自动回归、制造解、流动、MPI、性能 | 全套验收产物可复现，无手工修网格步骤，所有失败公开 |

基础几何基准使用工程内可重建的解析对象，不能只靠外部大 STEP；另外必须有至少一个非本程序导出的、来源和单位可追踪的 STEP 模型，防止只做自导出自导入测试。复杂模型可先用合成船体＋艉部有限厚尾缘，再加入公开且允许使用的船体/水翼 STEP，锁定文件hash。

建议数量目标：指定外流基准中，距物体及边界层足够远的背景核心保持严格六面体；全域严格hex数占比起始目标>=70%，单独报告其体积占比。此数只适用于锁定基准，不能对任意几何保证；若长期达不到应重新评估贴体/过渡路线，不许靠误标单元类型达标。

### 8.1 数值验收必须分层

1. 格式/几何：真实 `readMeshFile()` + `validate()`；核对cell数、体积、patch面积及形心。
2. 守恒/一致性：常数场、线性场重构、共享面通量等值反号、粗细界面无重复装配。先使用解析box/加密box隔离CAD误差。
3. 制造解：在3D box上选解析Poisson解，三档网格（尺寸比2），体积加权L2误差随细化下降，后两档观测阶>=1.5；局部加密族另做尺寸约束一致的误差研究。线性容差收紧10倍后误差变化<1%，证明不是迭代误差主导。失败需定位是网格还是求解离散，不直接篡改求解器通过。
4. 带BL层流：平行板Poiseuille，细网格速度体积L2相对误差<=1%，归一化进出口流量不平衡<=1e-6，三网格误差递减。明确压力梯度、黏度、边界及收敛条件。
5. 串行与MPI 1/2/4 ranks：比较同一全局cell上的场及积分量；基础双精度基准归一化差异目标<=1e-7，记录线性/外迭代容差并做敏感性检查。不声称位一致或扩展性已经证实。
6. 工程冒烟：船体/水翼只需在已支持的物理模型下自然完成受限步数，字段有限、无线性失败、满足约定残差及流量指标。该结果不能称为船舶阻力/自由面物理验证；更不能因设计了水线加密就推断BabelSim具备VOF求解能力。

若现有 BabelSim 公共接口不支持某一制造解测试条件，M0记录缺口，在 BabelMesh integration 目录编写独立测试驱动并链接公共API；需要改求解器时作为独立任务报告，不在本任务中偷偷放宽求解器验证。

### 8.2 自动化与性能

100k/500k/1M cells 阶梯测试，记录CPU/内存/编译器/版本/线程数、每阶段墙钟、峰值RSS、最终N/F/连接数、失败重试和导出成本。均匀核心与复杂CAD/BL分别测，避免以立方体性能冒充船体性能。

建议在固定8物理核/32GiB机器上把100万单元标准基准的初始工程预算设为30分钟、峰值RSS<=16GiB；这些是待实际测量检验的预算，不是预估成功承诺。M0锁定机器与基准，M8不许临时换硬件获得“通过”。先优化BVH查询、稀疏树和内存；CAD线程安全以具体API验证，不能任意并发修改TopoDS对象。

## 9. 主要风险与处理

| 风险 | 处理及停止条件 |
|---|---|
| STEP缝隙/单位错误/薄特征 | 相对与绝对容差、单位报告、修复前后体积；超容差失败 |
| 粗细界面T-junction | 共享面overlay和边环传播；逐cell边配对为硬门禁 |
| 贴体塌陷/跨薄壁误吸附 | 特征关联、候选限制、步长搜索、局部细化重试 |
| 极凹单元reader方向不兼容 | 优先凸单元、局部拆分/加密；不得修改reader绕过 |
| BL碰撞或覆盖不足 | 核心腾空间、碰撞预检、终止模板、strict失败 |
| 层顶与核心尺度突变 | 法向尺度匹配、局部核心加密、质量回滚 |
| 几何误差被“节点在CAD上”掩盖 | 检查面内双向采样、特征/通道完整性、体积偏差 |
| 自动化靠手工修复 | 从干净目录一条命令重放；所有参数变化归工程文件 |

最大技术风险是 M4 特征贴体和 M6 复杂边界层，不是 STEP 读入或均匀hex划分。应先做这些阶段的代表性基准并设置退出条件，再扩展更多CAD命令或GUI。实际工期要由 M3/M5 试点数据估算；不将工业网格器完整能力当作一次短开发即可完成的功能。

配套可执行任务文本见 [BabelMesh-开发Prompts.md](/home/midway/BabelMesh/docs/BabelMesh-开发Prompts.md)。所有命令和测试目标中未存在者均是后续开发要求，本轮没有执行开发。
