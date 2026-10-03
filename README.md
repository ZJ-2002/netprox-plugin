# netprox

Guney-Barabási 网络邻近（主计划 §18.2 CMap 候选的机制核验层）。
两组基因集在 PPI/代谢网络上的拓扑接近度 + 度匹配随机零模型。

## 节点

| kind | 输入 | 输出 |
|---|---|---|
| `network_proximity` | 边表 TSV（两列）+ 集合表 TSV（set_id/role/gene，role ∈ S,T） | `proximity.tsv`（set_s/set_t/measure/d_observed/d_rand_mean/d_rand_sd/z/p_lower/p_upper/n_s/n_t/n_random_used）、log |

参数：`measure`（closest|separation，默认 closest）、`n_random`
（默认 1000）、`bin_width`（默认 10）、`seed`。

## 布局

```text
netprox/
├── manifest.toml
├── scripts/network_proximity.py
├── Dockerfile
└── README.md
```

镜像：`localhost/autonomics/biotools-py@sha256:21e14c1582a9d4c261c0a23a11864293313ee2b49b9b7d85bcabc89210b48c7e`。

## 口径（Guney 2016 语义）

- **度匹配随机化**：从同度 bin（`degree // bin_width`）抽随机集；
  bin 不足时**聚合相邻 bin 至够用**，欠匹配计数每对集合只报一次
  （不随随机次数刷屏）。经验 P=(b+1)/(B+1) 双侧。
- `closest`：S 中每基因到 T 最近距离的均值；`separation`：
  d_AB − (d_AA + d_BB)/2，集合内距离**排除自身**（否则 d_AA≡0）。
- 只取最大连通分量；不在网络中的基因剔除并记 log（n_s/n_t 列报
  实际入网数）。
- BFS 按目标集合批量发起（每对 |T| 次 BFS，非每基因每目标两次）。

## 解释边界

网络邻近≠药效，≠直接互作；PPI 覆盖偏倚（well-studied 基因度更高）
本身是混杂。与 cmap 家族配合时：连通性筛候选，本节点只做机制
拓扑核验。

## 冒烟记录

稠密模块拓扑：近距对 d=1.75（经 hub ≤2 如设计）、远距对 z=+39.8
（P 下限 0.004975）；closest 与 separation 排序一致；度匹配随机
均值/SD 合理。
