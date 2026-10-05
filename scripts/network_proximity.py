# 中文注释：Guney 网络邻近（度数匹配随机化）。
# 输入构念：物理互作边表（同一命名空间）+ 种子集合表（set_id/role/gene）。
# 规则来源：Guney & Barabasi 2016（closest/separation，度数分箱匹配随机
# 对照）；方案 §17/§18.2——经验 P=(b+1)/(B+1)，双侧各报，置换数与种子
# 执行前锁死。统计口径：仅最大连通分量内的基因参与；分箱按度数除以
# bin_width 取整；箱内候选不足时记 undermatched 并如实输出（不降条件、
# 不聚合相邻箱）。
# v2 契约修正（审计 F12）：
# - 表头：has_header 参数（auto|yes|no，默认 auto）可显式覆盖。v1 的
#   header=None 会把 string_ppi 输出的表头当成一条假边（gene_a—gene_b）。
# - n_random 默认 1000 → 10000（对齐方案 §18.2 起始值；变更记录在 README）。
# - 输出列契约：实际输出列与 manifest 声明对齐为 n_random_used（实际
#   参与统计的随机次数），并新增 mc_se = d_rand_sd/sqrt(n_random_used)
#   （注意：这是随机距离均值的 MC 标准误，不是经验 P 的 MC 精度）。
#   v1 manifest 声明的 n_in_component 从未输出过。
# v3 契约修正（fix-review R05/R06，2026-10-05）：
# - 表头嗅探不再用"非数字且不在其余行"启发式——它会把只出现一次的真边
#   （如单边文件 "A\tB"）当表头删掉。auto 改为已知表头记号白名单：首行
#   两端点均为白名单记号（大小写不敏感）才判表头；任意基因名不承诺可
#   无歧义嗅探，非常规表头用 has_header=yes 显式声明。
# - separation 对齐官方定义（emreg00/toolbox get_separation, jorg-closest）：
#   d12 = 两个方向 nearest 距离均值的平均，separation = d12-(d_AA+d_BB)/2，
#   对称。单基因集合的集合内距离按 toolbox 以 0 代入（不再输出 NA）。
# - 失败保行：组合不可估计（集合经分量过滤后为空等）不再静默 skip 或
#   整体退出 2，输出 NA 行并 log 记因；仅当不存在任何 (S,T) 组合才退出 2。
# 失败处理：列缺失/无组合 → 退出 2。
# 解释边界：邻近=先验相容证据，不是药效（§18.2）。
import os
import sys
from collections import defaultdict

import networkx as nx
import numpy as np
import pandas as pd

DEFAULT_N_RANDOM = 10000  # v2：方案起始值；v1 为 1000
MEASURE = os.environ.get("NP_MEASURE", "closest").lower()
N_RANDOM = int(os.environ.get("NP_N_RANDOM", str(DEFAULT_N_RANDOM)))
BIN_WIDTH = int(os.environ.get("NP_BIN_WIDTH", "10"))
SEED = int(os.environ.get("NP_SEED", "42"))
HAS_HEADER = os.environ.get("NP_HAS_HEADER", "auto").lower()

OUTPUT_COLUMNS = [
    "set_s", "set_t", "measure", "d_observed", "d_rand_mean", "d_rand_sd",
    "mc_se", "z", "p_lower", "p_upper", "n_s", "n_t", "n_random_used",
]


def fail(message):
    print(message, file=sys.stderr)
    sys.exit(2)


def read_edge_rows(path):
    """读边表原始行：跳过 # 注释与空行，返回 (a, b) 端点对列表（第 3 列起忽略）。"""
    rows = []
    malformed = 0
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 2 or not fields[0].strip() or not fields[1].strip():
                malformed += 1
                continue
            rows.append((fields[0].strip(), fields[1].strip()))
    if not rows:
        fail(f"edge table has no parsable rows: {path}")
    if malformed:
        print(f"WARN skipped {malformed} malformed edge rows (<2 tab-separated columns)", file=sys.stderr)
    return rows


# 已知表头记号（fix-review R05）：本生态边表表头的实际用词，**小写精确
# 匹配**——大小写不敏感会把大写基因符号（如 GENE1/GENE2）误判成表头
# （gene1/gene2 在白名单内）。生态内表头均为小写 snake_case；大写或
# 非常规表头不承诺识别，需 has_header=yes 显式声明。单字母记号
# （a/b/u/v）也故意不收——真实基因符号/别名存在单字母形式。
HEADER_TOKENS = frozenset({
    "source", "target", "src", "dst", "from", "to",
    "node1", "node2", "gene1", "gene2", "gene_a", "gene_b",
    "protein1", "protein2", "interactor_a", "interactor_b",
    "partner1", "partner2", "id_a", "id_b", "symbol_a", "symbol_b",
    "geneid1", "geneid2", "preferred_name_a", "preferred_name_b",
})


def looks_like_header(rows):
    """首行两端点均为已知表头记号（小写精确匹配）→ 表头。

    不嗅探任意基因名：旧启发式（非数字且不出现在其余行端点）会把只出现
    一次的真边当表头删掉（单边文件 "A\tB" → 空图；两条不相连的符号边
    → 删第一条）。白名单外的首行一律当数据行。
    """
    first_a, first_b = rows[0]
    return first_a in HEADER_TOKENS and first_b in HEADER_TOKENS


def parse_edges(path, has_header="auto"):
    """返回 (边列表, 是否当作表头丢弃, 嗅探判定)。has_header: auto|yes|no。"""
    rows = read_edge_rows(path)
    if has_header == "yes":
        return rows[1:], True, "explicit"
    if has_header == "no":
        return rows, False, "explicit"
    if has_header != "auto":
        fail(f"has_header must be auto|yes|no, got {has_header}")
    sniffed = looks_like_header(rows)
    return (rows[1:] if sniffed else rows), sniffed, "sniffed"


def build_graph(edges):
    """返回 (整图, 最大连通分量节点集)。分量裁剪留给调用方（log 需要整图统计）。"""
    graph = nx.Graph()
    graph.add_edges_from(edges)
    component = max(nx.connected_components(graph), key=len)
    return graph, component


def closest_distance(graph, from_set, to_set, exclude_self=False):
    # 对 to_set 每个节点做 BFS，取 from_set 各点到 to_set 的最小距离均值。
    # 集合内距离（d_AA/d_BB）须排除自身，否则每点到自己的 0 会压平 d_AA；
    # 集合仅 1 个元素时排除自身后无可达目标 → None（separation 不可估计）。
    to_set = set(to_set)
    best = {node: None for node in from_set}
    for target in to_set:
        lengths = nx.single_source_shortest_path_length(graph, target)
        for node in best:
            if exclude_self and node == target:
                continue
            if node in lengths and (best[node] is None or lengths[node] < best[node]):
                best[node] = lengths[node]
    reachable = [value for value in best.values() if value is not None]
    if not reachable:
        return None
    return float(np.mean(reachable))


def proximity(graph, measure, s_members, t_members):
    """返回 (观测值, 不可估计原因)。原因 None = 正常估计。

    closest 保持 Guney 的方向性定义（S 中每基因到 T 最近距离的均值，
    S/T 互换会变——这是 closest 本身的性质，不是缺陷）。
    separation 对齐 emreg00/toolbox get_separation（jorg-closest）：
    d12 = (mean_{s∈S} min_{t∈T} d + mean_{t∈T} min_{s∈S} d) / 2，对称；
    separation = d12 − (d_SS + d_TT)/2，d_XX 为集合内 nearest 距离均值
    （排除自身）。单基因集合的 d_XX 按 toolbox 以 0 代入（values=[0]）。
    """
    d_st = closest_distance(graph, s_members, t_members)
    d_ts = closest_distance(graph, t_members, s_members)
    if d_st is None or d_ts is None:
        return None, "unreachable"  # 两集合在最大分量内本应连通，防御分支
    if measure == "closest":
        return d_st, None
    d_ss = closest_distance(graph, s_members, s_members, exclude_self=True)
    d_tt = closest_distance(graph, t_members, t_members, exclude_self=True)
    if d_ss is None:
        d_ss = 0.0  # 单基因集合：toolbox values=[0] 代入
    if d_tt is None:
        d_tt = 0.0
    d12 = (d_st + d_ts) / 2.0
    return d12 - 0.5 * (d_ss + d_tt), None


def degree_bins(graph, bin_width):
    bin_of = {node: int(degree // bin_width) for node, degree in graph.degree()}
    bins = defaultdict(list)
    for node, bucket in bin_of.items():
        bins[bucket].append(node)
    return bin_of, bins


def sample_matched(rng, bins, bin_of, original_members):
    # 返回 (抽样, 欠匹配箱数)。欠匹配只统计不刷 log，由调用方每对汇总一次。
    pool_by_bin = defaultdict(list)
    for node in original_members:
        pool_by_bin[bin_of[node]].append(node)
    chosen = []
    undermatched = 0
    for bucket, members in pool_by_bin.items():
        candidates = np.array(bins.get(bucket, []))
        take = min(len(members), len(candidates))
        picked = rng.choice(candidates, size=take, replace=False)
        chosen.extend(picked.tolist())
        if take < len(members):
            undermatched += 1
    return chosen, undermatched


def main():
    if MEASURE not in ("closest", "separation"):
        fail(f"measure must be closest|separation, got {MEASURE}")
    edges_path = os.environ["AUTONOMICS_INPUT0"]
    sets_path = os.environ["AUTONOMICS_INPUT1"]
    out_tsv = os.environ["AUTONOMICS_OUTPUT0"]
    out_log = os.environ["AUTONOMICS_OUTPUT1"]

    edge_rows, dropped_header, header_mode = parse_edges(edges_path, HAS_HEADER)
    full_graph, component = build_graph(edge_rows)
    graph = full_graph.subgraph(component).copy()
    bin_of, bins = degree_bins(graph, BIN_WIDTH)

    sets_table = pd.read_csv(sets_path, sep="\t")
    for column in ("set_id", "role", "gene"):
        if column not in sets_table.columns:
            fail(f"sets TSV missing column: {column}")
    sets = defaultdict(set)
    for _, row in sets_table.iterrows():
        sets[(str(row["set_id"]), str(row["role"]).upper())].add(str(row["gene"]))

    s_keys = sorted(key for key in sets if key[1] == "S")
    t_keys = sorted(key for key in sets if key[1] == "T")
    if not s_keys or not t_keys:
        fail("sets TSV needs at least one S-role set and one T-role set")

    rng = np.random.default_rng(SEED)
    log_lines = [
        # 整图统计（含分量外节点/边）+ 最大分量大小——表头假边之类的图污染在
        # nodes/edges 里可见，component 行说明实际参与统计的子图。
        f"nodes={full_graph.number_of_nodes()}\tedges={full_graph.number_of_edges()}\tcomponent={len(component)}",
        "header\tmode={}\tdropped={}  # auto sniff (v3): 首行两端点均为已知表头记号白名单成员才判表头；白名单外一律当数据行".format(
            header_mode, dropped_header
        ),
        "params\tmeasure={}\tn_random={}\tbin_width={}\tseed={}  # n_random default 10000 since v2 (was 1000); mc_se = MC SE of the null-mean distance, NOT of the empirical p".format(
            MEASURE, N_RANDOM, BIN_WIDTH, SEED
        ),
    ]

    def members_in_component(key):
        return sorted(sets[key] & component)

    records = []
    n_estimable = 0

    def na_row(s_key, t_key, s_members, t_members, reason):
        # 失败保行（fix-review R06）：不可估计的组合保留 NA 行（n_s/n_t 如实
        # 报分量内入网数），log 记因；不再静默 skip。
        nan = float("nan")
        records.append((s_key[0], t_key[0], MEASURE, nan, nan, nan,
                        nan, nan, nan, nan, len(s_members), len(t_members), nan))
        log_lines.append(f"NA {s_key[0]} vs {t_key[0]}: {reason}")

    for s_key in s_keys:
        for t_key in t_keys:
            s_members = members_in_component(s_key)
            t_members = members_in_component(t_key)
            if not s_members or not t_members:
                na_row(s_key, t_key, s_members, t_members,
                       "not estimable - a set is empty after largest-component filter "
                       "(no member genes in the largest component)")
                continue
            observed, reason = proximity(graph, MEASURE, s_members, t_members)
            if observed is None:
                na_row(s_key, t_key, s_members, t_members,
                       f"not estimable - {reason}")
                continue
            n_estimable += 1
            null = []
            undermatched_total = 0
            for _ in range(N_RANDOM):
                s_rand, under_s = sample_matched(rng, bins, bin_of, s_members)
                t_rand, under_t = sample_matched(rng, bins, bin_of, t_members)
                undermatched_total += under_s + under_t
                value, reason_rand = proximity(graph, MEASURE, s_rand, t_rand)
                if value is not None:
                    null.append(value)
            if undermatched_total > 0:
                log_lines.append(
                    f"WARN {s_key[0]} vs {t_key[0]}: degree bins undermatched "
                    f"{undermatched_total} times across {N_RANDOM} randomizations; "
                    "random sets are smaller than requested - treat as a deviation from "
                    "the frozen null model, not merely reduced precision (plan §17)"
                )
            null_array = np.asarray(null, dtype=float)
            mean = float(null_array.mean())
            sd = float(null_array.std(ddof=1)) if len(null_array) > 1 else 0.0
            z = (observed - mean) / sd if sd > 0 else 0.0
            p_lower = (1 + int((null_array <= observed).sum())) / (len(null_array) + 1)
            p_upper = (1 + int((null_array >= observed).sum())) / (len(null_array) + 1)
            n_used = len(null_array)
            mc_se = sd / np.sqrt(n_used) if n_used > 0 else float("nan")
            records.append((s_key[0], t_key[0], MEASURE, observed, mean, sd,
                            float(mc_se), z, p_lower, p_upper,
                            len(s_members), len(t_members), n_used))
            log_lines.append(
                f"done {s_key[0]} vs {t_key[0]}: z={z:.3f} p_lower={p_lower:.4f} B_used={n_used}"
            )

    # 失败保行（fix-review R06）：先写证据文件，再按可估性定退出码——全部
    # 组合不可估计 → 仍退 2（节点失败信号，证据已落盘）；部分可估 → 退 0，
    # NA 行随表保留。
    pd.DataFrame(records, columns=OUTPUT_COLUMNS).to_csv(out_tsv, sep="\t", index=False)
    with open(out_log, "w", encoding="utf-8") as handle:
        handle.write("\n".join(log_lines) + "\n")
    if n_estimable == 0:
        fail("no (S,T) combination was estimable on the largest component "
             "(all rows NA; evidence written to outputs)")


if __name__ == "__main__":
    main()
