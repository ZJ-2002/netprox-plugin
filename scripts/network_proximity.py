# 中文注释：Guney 网络邻近（度数匹配随机化）。
# 输入构念：物理互作边表（同一命名空间）+ 种子集合表（set_id/role/gene）。
# 规则来源：Guney & Barabasi 2016（closest/separation，度数分箱匹配随机
# 对照）；方案 §17/§18.2——经验 P=(b+1)/(B+1)，双侧各报，置换数与种子
# 执行前锁死。统计口径：仅最大连通分量内的基因参与；分箱按度数除以
# bin_width 取整；箱内候选不足时记 undermatched 并如实输出（不降条件）。
# 失败处理：列缺失/集合全在分量外/无组合 → 退出 2。
# 解释边界：邻近=先验相容证据，不是药效（§18.2）。
import os
import sys
from collections import defaultdict

import networkx as nx
import numpy as np
import pandas as pd

MEASURE = os.environ.get("NP_MEASURE", "closest").lower()
N_RANDOM = int(os.environ.get("NP_N_RANDOM", "1000"))
BIN_WIDTH = int(os.environ.get("NP_BIN_WIDTH", "10"))
SEED = int(os.environ.get("NP_SEED", "42"))
EDGES = os.environ["AUTONOMICS_INPUT0"]
SETS = os.environ["AUTONOMICS_INPUT1"]
OUT_TSV = os.environ["AUTONOMICS_OUTPUT0"]
OUT_LOG = os.environ["AUTONOMICS_OUTPUT1"]


def fail(message):
    print(message, file=sys.stderr)
    sys.exit(2)


if MEASURE not in ("closest", "separation"):
    fail(f"measure must be closest|separation, got {MEASURE}")

edges = pd.read_csv(EDGES, sep="\t", header=None, usecols=[0, 1],
                    names=["a", "b"], dtype=str, comment="#")
graph = nx.Graph()
graph.add_edges_from(zip(edges["a"], edges["b"]))
component = max(nx.connected_components(graph), key=len)
graph = graph.subgraph(component).copy()
degrees = dict(graph.degree())
bin_of = {node: int(degree // BIN_WIDTH) for node, degree in degrees.items()}
bins = defaultdict(list)
for node, bucket in bin_of.items():
    bins[bucket].append(node)

sets_table = pd.read_csv(SETS, sep="\t")
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
    f"nodes={graph.number_of_nodes()}\tedges={graph.number_of_edges()}\tcomponent={len(component)}",
    f"measure={MEASURE}\tn_random={N_RANDOM}\tbin_width={BIN_WIDTH}\tseed={SEED}",
]


def members_in_component(key):
    return sorted(sets[key] & component)


def closest_distance(from_set, to_set, exclude_self=False):
    # 对 to_set 每个节点做 BFS，取 from_set 各点到 to_set 的最小距离均值。
    # 集合内距离（d_AA/d_BB）须排除自身，否则每点到自己的 0 会压平 d_AA。
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


def proximity(s_members, t_members):
    d_st = closest_distance(s_members, t_members)
    if d_st is None:
        return None
    if MEASURE == "closest":
        return d_st
    d_ss = closest_distance(s_members, s_members, exclude_self=True)
    d_tt = closest_distance(t_members, t_members, exclude_self=True)
    return d_st - 0.5 * (d_ss + d_tt)


def sample_matched(original_members):
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


records = []
for s_key in s_keys:
    for t_key in t_keys:
        s_members = members_in_component(s_key)
        t_members = members_in_component(t_key)
        if not s_members or not t_members:
            log_lines.append(f"skipped {s_key[0]} vs {t_key[0]}: empty after component filter")
            continue
        observed = proximity(s_members, t_members)
        if observed is None:
            log_lines.append(f"skipped {s_key[0]} vs {t_key[0]}: unreachable")
            continue
        null = []
        undermatched_total = 0
        for _ in range(N_RANDOM):
            s_rand, under_s = sample_matched(s_members)
            t_rand, under_t = sample_matched(t_members)
            undermatched_total += under_s + under_t
            value = proximity(s_rand, t_rand)
            if value is not None:
                null.append(value)
        if undermatched_total > 0:
            log_lines.append(
                f"WARN {s_key[0]} vs {t_key[0]}: degree bins undermatched "
                f"{undermatched_total} times across {N_RANDOM} randomizations; "
                "random sets are smaller than requested - report as reduced precision (plan §17)"
            )
        null_array = np.asarray(null, dtype=float)
        mean = float(null_array.mean())
        sd = float(null_array.std(ddof=1)) if len(null_array) > 1 else 0.0
        z = (observed - mean) / sd if sd > 0 else 0.0
        p_lower = (1 + int((null_array <= observed).sum())) / (len(null_array) + 1)
        p_upper = (1 + int((null_array >= observed).sum())) / (len(null_array) + 1)
        records.append((s_key[0], t_key[0], MEASURE, observed, mean, sd, z,
                        p_lower, p_upper, len(s_members), len(t_members),
                        len(null_array)))
        log_lines.append(f"done {s_key[0]} vs {t_key[0]}: z={z:.3f} p_lower={p_lower:.4f}")

if not records:
    fail("no (S,T) combination was estimable on the largest component")

pd.DataFrame(records, columns=[
    "set_s", "set_t", "measure", "d_observed", "d_rand_mean", "d_rand_sd",
    "z", "p_lower", "p_upper", "n_s", "n_t", "n_random_used",
]).to_csv(OUT_TSV, sep="\t", index=False)
with open(OUT_LOG, "w", encoding="utf-8") as handle:
    handle.write("\n".join(log_lines) + "\n")
