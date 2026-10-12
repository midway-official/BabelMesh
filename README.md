# BabelMesh

BabelMesh 提供 C++17 网格核心、命令行与 Qt 6 CAD/网格桌面界面。新工程使用 `schema_version: 2`：**先生成正交笛卡尔背景，再挖除非流体与预留区域，最后构造边界层和贴体混合单元**。保留 OCCT CAD、STEP、旧版笛卡尔/Octree 和历史验收资料；新主线不再整体吸附背景点。

| 生成规则 | 笛卡尔区域 | 壁面及贴体区域 | 最终输出 |
|---|---|---|---|
| 二维 | 四边形挤出为六面体；加密接口可为多边形柱 | 三角形挤出为三棱柱 | 一层体网格，前后面 symmetry；禁止四面体/四棱锥 |
| 三维 | 正交六面体；粗细接口预留为混合过渡区 | 四面体；六面体接口按需使用四棱锥 | 三维混合体网格 |

共享面和可变单元连接是真实数据。输出采用求解器当前的 `BABELSIM_MESH 3` 面格式；测试直接编译相邻 `../BabelSim` 的 reader 核对拓扑与几何。

## 构建和运行

Ubuntu 22.04 依赖安装示例：

```bash
sudo apt install make cmake g++ python3 python3-jsonschema \
  nlohmann-json3-dev libeigen3-dev gmsh libgmsh-dev libcgns-dev \
  qt6-base-dev libqt6opengl6-dev libqt6svg6 \
  libocct-foundation-dev libocct-modeling-data-dev \
  libocct-modeling-algorithms-dev libocct-ocaf-dev \
  libocct-visualization-dev libocct-data-exchange-dev

cmake -S . -B build/qt -DCMAKE_BUILD_TYPE=Release
cmake --build build/qt -j4
build/qt/bin/babelmesh-gui examples/cylinder2d.json
```

当前工作区已构建上述程序。已有 OCCT/Qt 基础依赖的 Jammy 环境可用 [开发依赖缓存脚本](tools/bootstrap_hybrid_dependencies.sh) 补充本地 Gmsh/CGNS 开发包，无需管理员权限。仅构建 CLI 可设置 `-DBABELMESH_BUILD_GUI=OFF`；reader 路径可用 `-DBABELSIM_ROOT=/absolute/path/BabelSim` 指定。

```bash
build/qt/bin/babelmesh generate examples/cylinder2d.json --out outputs/cylinder2d
build/qt/bin/babelmesh generate examples/sphere3d.json --out outputs/sphere3d
build/qt/bin/babelmesh check outputs/cylinder2d/mesh.json
build/qt/bin/babelmesh optimize outputs/cylinder2d/mesh.json --out outputs/optimized
build/qt/bin/babelmesh export-cad examples/sphere3d.json --entity body --out outputs/body.step
```

默认保护已有结果；明确需要覆盖时使用 `--overwrite` 或工程中的 `output.overwrite: true`。生成和检查完成后才发布文件；严格模式失败不发布新网格。报告分别表示拓扑通过、面质量通过和诊断结果。**生成成功、reader 通过与完整 CFD 质量验收是不同结论。**

## CAD 和界面

界面提供 STEP 导入、基本体、多边形草图/拉伸、布尔、变换、旋转成形和放样，操作树和工程 JSON、撤销/重做、拾取面并分配物理组、旋转/平移/缩放、轴视图、真实网格颜色和裁剪、参数设置、后台任务与取消。复杂操作通过参数对话框编辑。“网格”菜单调用转换和优化，“文件”菜单导出流体域 STEP。

预览来自实际核心输出。任务使用工程快照；生成期间继续修改工程，结果会标记为启动时的参数。工程长度单位与米制输出分开，CAD 贴合容差单独用米指定。

视图提供正交/透视切换、磨砂半透明材质、不透明度和黑色实线单元边。三维默认使用真实平面剖切，生成共面的单元截面；可选择 X/Y/Z 剖切或整单元筛选。壁面用中性色区分，曲面着色保留尖角。显示间隙只影响完整单元视图，平面截面不收缩，不修改求解器网格。

## 参数、质量与交换

背景可指定 `cell_counts` 或 `cell_size`。完整覆盖的外盒边界直接保留笛卡尔单元，避免无必要的外圈贴体带。局部 box/sphere 尺寸源通过 `apply_to` 控制 background、surface、transition 或 all。二维采用 quadtree，三维采用 octree；粗细邻接采用 2:1 约束。边界层指定目标物理组、首层高度、层数和增长率；报告核对实际法向层高与覆盖。

背景与计算域一旦指定，生成、贴体、层带和优化均不得越界。二维截面加减半厚度也必须落在背景及计算域中；CAD 域超出背景会明确拒绝，不自动扩大背景。最终报告检查全部节点；优化的每个候选使用同一硬约束。通过本工具导出的 Fluent/CGNS 保留锁定边界元数据。

自动优化固定背景、物理边界和边界层节点，二维成对移动前后节点；三维默认尝试少量 3→2、4→4 过渡四面体换边，再按新邻接平滑。可选合并少量过渡四面体为真实凸多面体。检查正体积与凸壳并回滚较差方案，**增加单元数为 0**。表面离散为满足 CAD 容差而发生的有限加密属于生成阶段，受尺寸下限与资源预算控制，单独报告。论文依据和实现选择见 [高质量混合网格论文与实现](docs/高质量混合网格论文与实现.md)。

继续优化优先使用 `mesh.json` 检查点。原生 `.bsm` 没有生成区域标签；重新读入时可恢复严格成对的二维挤出规则并识别轴对齐背景，但无法恢复边界层保护信息。因此缺失区域标签的导入文件仅检查和报告，保留全部节点和单元；已保存完整标签的检查点可执行节点优化。

当前已实现闭壳、共享面、边方向、边界归属、正体积、凸壳、非正交角、内部偏斜、双向采样 CAD 距离、实际层高与层路径间隙检查。全局单元交叠扫描和历史 CFD profile 的所有指标尚未评估，报告保留 `full_cfd_profile_evaluated: false`。三维加密或边界层案例可能只通过拓扑而未通过面质量，应按报告调整尺寸/层参数；严格模式拒绝相关结果。

```bash
build/qt/bin/babelmesh convert outputs/cylinder2d/mesh.json --to fluent-ascii --out outputs/grid.msh
build/qt/bin/babelmesh convert outputs/grid.msh --to babelsim-v3 --out outputs/imported.bsm
build/qt/bin/babelmesh convert outputs/cylinder2d/mesh.json --to cgns --out outputs/grid.cgns
build/qt/bin/babelmesh convert outputs/grid.cgns --to checkpoint --out outputs/imported.json
```

双向交换支持 Fluent ASCII 三维面网格、CGNS 单一非结构体网格区、BabelSim v3 和检查点；VTU 支持导出。CGNS 支持标准混合单元和 NGON/NFACE 任意多面体（`cgns-poly`）。转换坐标、连接与边界组，不转换物理模型、数值设置或场数据。Fluent 面法向按 c0 约定转换，CGNS 读取长度单位并归一化为米。

STAR-CCM+ 和 FINE/Marine 需从原软件导出受支持的交换文件；本次未实现 `.ccm/.sim/.hex/.igg` 原生数据库，也未在商业软件内进行导入验收。FINE/Marine 的整网格交换选项需在其导入/导出界面确认，不能把 CGNS 曲面接口当成任意体网格接口。Fluent binary/HDF5、原生平面二维、多区/overset/周期连接输入会给出明确的不支持诊断。

参数的单位、约束、实例和格式边界见 [参数与格式说明](docs/混合网格参数与格式说明.md)，架构决定见 [重构方案](docs/混合网格重构方案.md)。

## 验证和历史路径

```bash
make test-hybrid
make MODE=release test -j4
```

新集成测试包括真实二维/三维 CAD、布尔和 STEP、三个截面轴、局部加密、边界层、格式往返、真实 BabelSim 几何核对、资源限制、覆盖保护与严格质量失败。可用时额外运行独立 CGNS cgnscheck。Qt 的 `--workflow-test` 验证实际编辑、保存、后台生成和渲染。

本轮复核目录为 `acceptance/evidence/hybrid-20261012`，详细结论见 [重构交付与验收说明](acceptance/HYBRID_REWRITE_20261012.md)；2026-10-10 的早期结果也保留。`schema_version: 1` 仍走保留的旧路径，其尺寸语义、M0–M4 门槛和命令见 [旧版说明](docs/旧版M0-M4使用说明.md)。旧曲面证据与新算法分开记录。依赖版本和许可证见 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。
