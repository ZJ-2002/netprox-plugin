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

- **表头处理（v2，审计 F12）**：边表首行自动嗅探——两端点列均非数字
  且不出现在其余任何行的端点集里 → 判为表头丢弃。v1 `header=None`
  会把 string_ppi 输出的 `gene_a	gene_b` 表头当成一条假边入图。
  `has_header=yes|no` 显式覆盖（no 时表头行按数据行处理，自负）。
- **度匹配随机化**：从同度 bin（`degree // bin_width`）无放回抽样；
  箱内候选不足时如实少抽并计 undermatched（每对集合汇总一次上报），
  **不降条件、不聚合相邻箱**（此处修正 v1 README 的"聚合相邻 bin"
  错误描述——代码从未聚合，聚合会放松度匹配）。经验 P=(b+1)/(B+1)
  双侧。
- `closest`：S 中每基因到 T 最近距离的均值；`separation`：
  d_AB − (d_AA + d_BB)/2，集合内距离**排除自身**（否则 d_AA≡0）。
- **单元素集合（v2）**：separation 某集合入网节点 <2 时，排除自身后
  集合内距离无定义 → 该 (S,T) 行输出 NA（全部统计列 NA，n_s/n_t 仍
  报实际入网数），不崩溃；log 记 "separation not estimable"。
- 只取最大连通分量；不在网络中的基因剔除并记 log（n_s/n_t 列报
  实际入网数）。
- **输出列契约（v2）**：`n_random_used` = 实际参与统计的随机化次数
  （v1 manifest 声明的 `n_in_component` 从未真实输出，属契约缺口，
  已改为声明实际列并由测试锁定）；`mc_se` = d_rand_sd/√n_random_used
  （MC 标准误，v2 新增）。log 逐对打印实际 B_used。
- BFS 按目标集合批量发起（每对 |T| 次 BFS，非每基因每目标两次）。

## 解释边界

网络邻近≠药效，≠直接互作；PPI 覆盖偏倚（well-studied 基因度更高）
本身是混杂。与 cmap 家族配合时：连通性筛候选，本节点只做机制
拓扑核验。**表达/可检测性匹配不在本节点**（plan_network_background
任务范围，勿在此扩）。

## 测试（离线，tests/test_network_proximity.py）

- 金标准：合成小图上用**纯 python BFS**（不经 networkx）独立算
  closest/separation，与脚本函数逐一比对（防实现两侧同错）。
- 表头两态：同一图有/无表头跑出的 proximity.tsv **逐字节一致**
  （假边不再入图）；`has_header=no` 显式覆盖时表头按数据行入图
  （log 节点数可证）。
- 单元素 separation → NA 行 + 退出 0 + log 记因；closest 不受影响。
- seed 确定性：同 seed 两次运行输出逐字节一致。
- 契约锁定：输出列 == manifest doc 声明列（tomllib 解析比对）；
  默认 n_random=10000；p ∈ [1/(B+1), 1]。

## 冒烟记录

稠密模块拓扑：近距对 d=1.75（经 hub ≤2 如设计）、远距对 z=+39.8
（P 下限 0.004975）；closest 与 separation 排序一致；度匹配随机
均值/SD 合理。
