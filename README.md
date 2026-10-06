# netprox

Guney-Barabási 网络邻近（主计划 §18.2 CMap 候选的机制核验层）。
两组基因集在 PPI/代谢网络上的拓扑接近度 + 度匹配随机零模型。

## 节点

| kind | 输入 | 输出 |
|---|---|---|
| `network_proximity` | 边表 TSV（两列）+ 集合表 TSV（set_id/role/gene，role ∈ S,T） | `proximity.tsv`（set_s/set_t/measure/d_observed/d_rand_mean/d_rand_sd/mc_se/z/p_lower/p_upper/n_s/n_t/n_random_used）、log |

参数：`measure`（closest|separation，默认 closest）、`n_random`
（默认 **10000**，v2 起对齐方案 §18.2 起始值——v1 默认 1000，变更于
2026-10-05 审计 F12）、`bin_width`（默认 10）、`seed`（默认 42）、
`has_header`（v2 新增：auto|yes|no，默认 auto）。

## 布局

```text
netprox/
├── manifest.toml
├── scripts/network_proximity.py
├── tests/test_network_proximity.py   # 离线：金标准 + 契约锁定
├── Dockerfile
└── README.md
```

镜像：`localhost/autonomics/biotools-py@sha256:21e14c1582a9d4c261c0a23a11864293313ee2b49b9b7d85bcabc89210b48c7e`。

## 口径（Guney 2016 语义）

- **表头处理（v3，fix-review R05）**：auto 模式改为**已知表头记号白名单**
  ——首行两端点均为白名单记号（`gene_a`/`gene_b`/`source`/`target`/
  `node1`/`protein1` 等，**小写精确匹配**）才判表头丢弃；白名单外一律当
  数据行。v2 的"非数字且不出现在其余行端点集"启发式会把只出现一次的
  真边当表头删掉（单边文件 `A	B` → 空图；两条不相连符号边 → 删
  第一条）。单字母记号（a/b/u/v）故意不收——真实基因符号存在单字母
  形式；大写形式（`SOURCE`/`GENE1`）也不识别——大写是真实基因符号的
  常态，大写或非常规表头请 `has_header=yes` 显式声明。v1 `header=None`
  会把 string_ppi 输出的表头当假边入图（v2 修复的原始缺陷）。
- **度匹配随机化**：从同度 bin（`degree // bin_width`）无放回抽样；
  箱内候选不足时如实少抽并计 undermatched（每对集合汇总一次上报），
  **不降条件、不聚合相邻箱**（此处修正 v1 README 的"聚合相邻 bin"
  错误描述——代码从未聚合，聚合会放松度匹配）。经验 P=(b+1)/(B+1)
  双侧。欠匹配抽出的随机集合比请求的小——按偏离冻结零模型对待并
  WARN，不只当精度损失。
- `closest`：S 中每基因到 T 最近距离的均值（方向性是 closest 的定义
  性质，S/T 互换会变）；`separation`（v4，2026-10-06 复审修正）对齐
  emreg00/toolbox `get_separation`（jorg-closest）**原文**：d12 = 两个
  方向**逐节点** nearest 距离合并成一个列表后的**整体均值**（toolbox 先
  `values.extend(...)` 再取一次 mean），separation = d12 − (d_AA +
  d_BB)/2，集合内距离**排除自身**（否则 d_AA≡0），**对称**。v3 曾错写成
  "两方向均值的等权平均"——|S|≠|T| 时两者不同（判别性反例：路径图
  A–B–C–D–E–F–G–H 上 S={A,B}/T={B,E,H}，v3 得 −0.25，toolbox 原文得
  0.00；v3 的反例 S={A,B}/T={B,H} 恰好 |S|=|T|，两口径同值 −1.75，
  不具判别力，v2 单向实现则得 −3.0/−0.5）。
- **单元素集合（v3）**：separation 的集合内距离按 toolbox 以 0 代入
  （`values=[0]`），单基因集合可估计（v2 曾输出 NA，官方语义可估）。
- **失败保行（v3）**：组合不可估计（某集合经最大分量过滤后为空等）
  保留 NA 行（统计列 NA，n_s/n_t 报实际入网数）+ log 记因；v2 曾静默
  skip 该组合。**全部组合**均不可估计 → 证据文件照写后退 2（节点失败
  信号不丢）。
- 只取最大连通分量；不在网络中的基因剔除并记 log（n_s/n_t 列报
  实际入网数）。
- **输出列契约（v2）**：`n_random_used` = 实际参与统计的随机化次数
  （v1 manifest 声明的 `n_in_component` 从未真实输出，属契约缺口，
  已改为声明实际列并由测试锁定）；`mc_se` = d_rand_sd/√n_random_used
  ——**随机距离均值的 MC 标准误，不是经验 P 的 MC 精度**（v3 措辞
  对齐）。log 逐对打印实际 B_used。
- BFS 按目标集合批量发起（每对 |T| 次 BFS，非每基因每目标两次）。

## 解释边界

网络邻近≠药效，≠直接互作；PPI 覆盖偏倚（well-studied 基因度更高）
本身是混杂。与 cmap 家族配合时：连通性筛候选，本节点只做机制
拓扑核验。**表达/可检测性匹配不在本节点**（plan_network_background
任务范围，勿在此扩）。

## 测试（离线，tests/test_network_proximity.py）

- 金标准：合成小图上用**纯 python BFS**（不经 networkx）独立算
  closest/separation（v4 金标镜像 toolbox 作者代码路径：两方向逐节点
  值 pooled + 单基因集合内距离 0；v3 金标照抄了插件公式，两侧同错
  一直绿），与脚本函数逐一比对，并用硬数值断言锁口径。
- 表头两态：同一图有/无表头跑出的 proximity.tsv **逐字节一致**
  （假边不再入图）；`has_header=no` 显式覆盖时表头按数据行入图
  （log 节点数可证）。
- 反例锁定：单边无表头文件 `A	B` 不再被嗅探删除；两条不相连
  符号边两条全保留；路径图 S={A,B}/T={B,H} separation 两向对称
  （= −1.75，v2 得 −3.0/−0.5）；**v4 判别性反例** S={A,B}/T={B,E,H}
  （|S|≠|T|）separation = 0.00（v3 错口径得 −0.25，此断言在 v3 代码
  上必挂）；单基因集合可估计（pooled 口径硬数值 2.75）；空集合（分量
  过滤后）→ NA 行保留不跳过。
- seed 确定性：同 seed 两次运行输出逐字节一致。
- 契约锁定：输出列 == manifest doc 声明列（tomllib 解析比对）；
  默认 n_random=10000；p ∈ [1/(B+1), 1]。

## 冒烟记录

（v3 代码记录，2026-10-05。）稠密模块拓扑：近距对 d=1.75（经 hub ≤2
如设计；此为 closest 值，v4 未动 closest）、远距对 z=+39.8（P 下限
0.004975）；closest 与 separation 排序一致；度匹配随机均值/SD 合理。
v4 只改 separation 在 |S|≠|T| 时的数值口径（向 toolbox 原文对齐），
seed 固定下 separation 的 z/p 需以 v4 重跑为准。
