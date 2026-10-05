# 离线测试（审计 F12）：金标准（纯 python BFS 独立实现）+ 表头两态 +
# 单元素 NA + seed 确定性 + 输出列契约锁定。不打网络、不依赖容器。
import importlib.util
import os
import subprocess
import sys
import tomllib
from collections import deque
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "network_proximity.py"
MANIFEST = REPO / "manifest.toml"
PYTHON = sys.executable


def load_module():
    spec = importlib.util.spec_from_file_location("network_proximity", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules["network_proximity"] = module
    spec.loader.exec_module(module)
    return module


# --- 独立金标准：纯 python BFS，不经 networkx/numpy（防两侧同错）。---

def bfs_distances(adjacency, source):
    distances = {source: 0}
    queue = deque([source])
    while queue:
        node = queue.popleft()
        for neighbor in adjacency[node]:
            if neighbor not in distances:
                distances[neighbor] = distances[node] + 1
                queue.append(neighbor)
    return distances


def all_distances(adjacency):
    return {node: bfs_distances(adjacency, node) for node in adjacency}


def closest_golden(dist, members, targets, exclude_self=False):
    values = []
    for node in members:
        candidates = [
            dist[node][t] for t in targets
            if t in dist[node] and not (exclude_self and t == node)
        ]
        if candidates:
            values.append(min(candidates))
    return sum(values) / len(values) if values else None


def separation_golden(dist, s_members, t_members):
    # v3（fix-review R06）：对齐 emreg00/toolbox get_separation（jorg-closest）
    # ——d12 取两个方向 nearest 距离均值的平均（对称），单基因集合的
    # 集合内距离以 0 代入（toolbox values=[0]）。v2 金标只有 S→T 单向。
    d_ab = closest_golden(dist, s_members, t_members)
    d_ba = closest_golden(dist, t_members, s_members)
    d_aa = closest_golden(dist, s_members, s_members, exclude_self=True)
    d_bb = closest_golden(dist, t_members, t_members, exclude_self=True)
    if d_ab is None or d_ba is None:
        return None
    if d_aa is None:
        d_aa = 0.0
    if d_bb is None:
        d_bb = 0.0
    return (d_ab + d_ba) / 2 - 0.5 * (d_aa + d_bb)


# 合成小图：hub 中心 + 两条分支（S 与 T 各挂一支）+ 一条孤立二点组件（应被剔除）。
GOLDEN_EDGES = [
    ("HUB", "S1"), ("HUB", "S2"), ("HUB", "T1"), ("HUB", "T2"),
    ("S1", "S2"), ("T1", "T3"), ("T2", "T3"), ("S1", "X1"), ("T1", "X2"),
    ("HUB", "FAR1"), ("FAR1", "FAR2"), ("FAR2", "FAR3"),  # 远端支
    ("ISO1", "ISO2"),  # 小组件：整组件剔除
]
S_SET = ["S1", "S2", "X1"]
T_SET = ["T1", "T2", "T3"]
T_SINGLETON = ["T3"]


def adjacency_of(edges):
    adjacency = {}
    for a, b in edges:
        adjacency.setdefault(a, set()).add(b)
        adjacency.setdefault(b, set()).add(a)
    return adjacency


def golden_component(adjacency):
    # 最大连通分量（独立实现：逐节点 BFS 分组取最大）。
    seen, best = set(), set()
    for node in adjacency:
        if node in seen:
            continue
        group = set(bfs_distances(adjacency, node))
        seen |= group
        if len(group) > len(best):
            best = group
    return best


GOLDEN_ADJ = adjacency_of(GOLDEN_EDGES)
GOLDEN_COMP = golden_component(GOLDEN_ADJ)
GOLDEN_DIST = all_distances({n: GOLDEN_ADJ[n] for n in GOLDEN_COMP})
S_IN = [g for g in S_SET if g in GOLDEN_COMP]
T_IN = [g for g in T_SET if g in GOLDEN_COMP]


def test_closest_matches_golden():
    mod = load_module()
    import networkx as nx

    graph = nx.Graph()
    graph.add_edges_from(GOLDEN_EDGES)
    component = max(nx.connected_components(graph), key=len)
    graph = graph.subgraph(component).copy()
    expected = closest_golden(GOLDEN_DIST, S_IN, T_IN)
    assert mod.closest_distance(graph, S_IN, T_IN) == expected


def test_separation_matches_golden_and_symmetric():
    mod = load_module()
    import networkx as nx

    graph = nx.Graph()
    graph.add_edges_from(GOLDEN_EDGES)
    component = max(nx.connected_components(graph), key=len)
    graph = graph.subgraph(component).copy()
    expected = separation_golden(GOLDEN_DIST, S_IN, T_IN)
    value, reason = mod.proximity(graph, "separation", S_IN, T_IN)
    assert value == expected and reason is None
    # 对称性（fix-review R06）：交换 S/T 集合，separation 不变（v2 单向实现
    # 会变）。closest 的方向性是其定义性质，不在此约束内。
    swapped, _ = mod.proximity(graph, "separation", T_IN, S_IN)
    assert swapped == value
    # 单元素集合（v3）：toolbox values=[0] 代入集合内距离 → 可估计，
    # 与独立金标一致（v2 输出 NA）。
    expected_singleton = separation_golden(GOLDEN_DIST, S_IN, T_SINGLETON)
    value_singleton, reason_singleton = mod.proximity(
        graph, "separation", S_IN, T_SINGLETON)
    assert reason_singleton is None
    assert value_singleton == expected_singleton


def test_separation_path_graph_counterexample_symmetric():
    # fix-review R06 原始反例：路径图 A–B–C–D–E–F–G–H，S={A,B}、T={B,H}。
    # v2 单向实现：separation(S,T)=−3.0、separation(T,S)=−0.5（非对称）。
    # 官方（toolbox jorg-closest）：d_ST=(1+0)/2=0.5、d_TS=(0+6)/2=3.0、
    # d_SS=1、d_TT=6 → (0.5+3.0)/2−(1+6)/2 = **−1.75**，两向一致。
    mod = load_module()
    import networkx as nx

    nodes = list("ABCDEFGH")
    graph = nx.Graph()
    graph.add_edges_from(zip(nodes, nodes[1:]))
    forward, _ = mod.proximity(graph, "separation", ["A", "B"], ["B", "H"])
    backward, _ = mod.proximity(graph, "separation", ["B", "H"], ["A", "B"])
    assert forward == -1.75 and backward == -1.75


def test_header_sniff_drops_gene_header_row(tmp_path):
    mod = load_module()
    headered = tmp_path / "edges.tsv"
    headered.write_text("gene_a\tgene_b\n" + EDGES_TSV, encoding="utf-8")
    edges, dropped, mode = mod.parse_edges(str(headered), "auto")
    assert dropped is True and mode == "sniffed"
    assert ("gene_a", "gene_b") not in edges
    # 无表头纯数据：不误删首行。
    plain = tmp_path / "plain.tsv"
    plain.write_text(EDGES_TSV, encoding="utf-8")
    edges, dropped, mode = mod.parse_edges(str(plain), "auto")
    assert dropped is False and mode == "sniffed" and len(edges) == len(GOLDEN_EDGES)


def test_header_sniff_keeps_single_symbolic_edge_review_r05(tmp_path):
    # fix-review R05 原始反例：真实无表头输入仅一条 "A\tB"。v2 启发式
    # （首行非数字且不出现在其余行端点）→ 误判表头 → 空边集。v3 白名单：
    # A/B 不是已知表头记号 → 保留为真边。
    mod = load_module()
    single = tmp_path / "single.tsv"
    single.write_text("A\tB\n", encoding="utf-8")
    edges, dropped, mode = mod.parse_edges(str(single), "auto")
    assert dropped is False and mode == "sniffed"
    assert edges == [("A", "B")]


def test_header_sniff_keeps_both_disconnected_symbolic_edges_review_r05(tmp_path):
    # fix-review R05 第二反例：两条互不相连的符号边——v2 会删第一条
    # （其端点不再出现于其余行）。v3 两条都保留。
    mod = load_module()
    two = tmp_path / "two.tsv"
    two.write_text("GENE1\tGENE2\nGENE3\tGENE4\n", encoding="utf-8")
    edges, dropped, mode = mod.parse_edges(str(two), "auto")
    assert dropped is False and mode == "sniffed"
    assert edges == [("GENE1", "GENE2"), ("GENE3", "GENE4")]


def test_header_sniff_whitelist_recognizes_known_headers(tmp_path):
    # 生态内真实表头（string_ppi 输出 gene_a/gene_b；通用 source/target）
    # 仍被识别——**小写精确匹配**。大写形式（SOURCE/Gene_A）与白名单词的
    # 大写基因符号（GENE1/GENE2 是 gene1/gene2 的大写）都不识别：
    # 大写或非常规表头需 has_header=yes 显式声明。
    mod = load_module()
    for header in ("gene_a\tgene_b", "source\ttarget", "node1\tnode2",
                   "protein1\tprotein2", "from\tto"):
        path = tmp_path / f"h{abs(hash(header))}.tsv"
        path.write_text(header + "\nA\tB\n", encoding="utf-8")
        edges, dropped, mode = mod.parse_edges(str(path), "auto")
        assert dropped is True and mode == "sniffed", header
        assert edges == [("A", "B")]
    for not_header in ("SOURCE\tTARGET", "Gene_A\tGene_B", "GENE1\tGENE2"):
        path = tmp_path / f"n{abs(hash(not_header))}.tsv"
        path.write_text(not_header + "\nA\tB\n", encoding="utf-8")
        edges, dropped, mode = mod.parse_edges(str(path), "auto")
        assert dropped is False and mode == "sniffed", not_header
        assert ("A", "B") in edges and (not_header.split("\t")[0],
                                        not_header.split("\t")[1]) in edges


# --- 全流程：子进程跑脚本（与引擎同构：env in / 文件 out）。---

EDGES_TSV = "\n".join("\t".join(row) for row in GOLDEN_EDGES) + "\n"
SETS_TSV = (
    "set_id\trole\tgene\n"
    + "".join(f"SSET\tS\t{g}\n" for g in S_SET)
    + "".join(f"TSET\tT\t{g}\n" for g in T_SET)
)


def run_node(tmp_path, edges_text, sets_text, extra_env=None, expect_ok=True):
    tmp_path.mkdir(parents=True, exist_ok=True)
    edges = tmp_path / "edges.tsv"
    sets = tmp_path / "sets.tsv"
    out_tsv = tmp_path / "proximity.tsv"
    out_log = tmp_path / "netprox.log"
    edges.write_text(edges_text, encoding="utf-8")
    sets.write_text(sets_text, encoding="utf-8")
    env = dict(os.environ)
    env.update({
        "AUTONOMICS_INPUT0": str(edges),
        "AUTONOMICS_INPUT1": str(sets),
        "AUTONOMICS_OUTPUT0": str(out_tsv),
        "AUTONOMICS_OUTPUT1": str(out_log),
        "NP_N_RANDOM": "100",
        "NP_BIN_WIDTH": "3",
    })
    env.update(extra_env or {})
    result = subprocess.run(
        [PYTHON, str(SCRIPT)], env=env, capture_output=True, text=True,
    )
    if expect_ok:
        assert result.returncode == 0, result.stderr
    return out_tsv, out_log, result


def read_table(path):
    # 注意不能用 .strip()：NA 行以 \t 结尾（pandas NaN→空串），strip 会吃掉
    # 行尾 tab 导致最后一列丢失。
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in lines[1:]]
    return header, rows


def test_headered_and_headerless_edge_tables_give_identical_output(tmp_path):
    headered = "gene_a\tgene_b\n" + EDGES_TSV
    plain = EDGES_TSV
    out1, _, _ = run_node(tmp_path / "a", headered, SETS_TSV)
    out2, _, _ = run_node(tmp_path / "b", plain, SETS_TSV)
    # v1 的 header=None 会把表头当假边（节点数+2、边+1）；v2 嗅探后两态一致。
    assert out1.read_bytes() == out2.read_bytes()
    # 与金标准比对：d_observed 逐位一致（度匹配抽样为随机项，另验分布界）。
    _, rows = read_table(out1)
    assert len(rows) == 1
    assert float(rows[0]["d_observed"]) == closest_golden(GOLDEN_DIST, S_IN, T_IN)
    assert rows[0]["n_s"] == str(len(S_IN)) and rows[0]["n_t"] == str(len(T_IN))


def test_has_header_override_semantics(tmp_path):
    headered = "gene_a\tgene_b\n" + EDGES_TSV
    out_auto, log_auto, _ = run_node(tmp_path / "auto", headered, SETS_TSV)
    out_no, log_no, _ = run_node(
        tmp_path / "no", headered, SETS_TSV, extra_env={"NP_HAS_HEADER": "no"},
    )
    # 大图上假边落在最大分量之外、数值恰好不变，但图被污染（log 整图统计
    # 自锚）：auto 丢弃表头 → 13 节点 13 边；显式 no → 假节点入图 → 15/14。
    assert "nodes=13\tedges=13\tcomponent=11" in log_auto.read_text()
    assert "nodes=15\tedges=14\tcomponent=11" in log_no.read_text()


def test_header_edge_can_hijack_largest_component_v1_failure_mode(tmp_path):
    # v1 真正会坏的场景：网络只有一条真边（A-B，两点分量）。表头假边
    # gene_a-gene_b 同为两点分量且先入图 → max(connected_components) 取到
    # 假分量 → A/B 全在分量外 → 整体不可估计退出 2。v2 auto 嗅探修复。
    tiny_edges = "gene_a\tgene_b\nA\tB\n"
    tiny_sets = "set_id\trole\tgene\nSSET\tS\tA\nTSET\tT\tB\n"
    out_tsv, _, result = run_node(tmp_path / "auto", tiny_edges, tiny_sets)
    assert result.returncode == 0, result.stderr
    _, rows = read_table(out_tsv)
    assert float(rows[0]["d_observed"]) == 1.0  # A→B 直连
    # 显式 has_header=no 恢复 v1 语义 → 复现失败模式（退出 2）。
    _, _, result_no = run_node(
        tmp_path / "no", tiny_edges, tiny_sets,
        extra_env={"NP_HAS_HEADER": "no"}, expect_ok=False,
    )
    assert result_no.returncode == 2


def test_separation_singleton_estimable_toolbox_zero(tmp_path):
    # v3（fix-review R06）：单基因集合按官方 toolbox values=[0] 代入集合内
    # 距离 → 可估计（v2 输出 NA 行）。数值与独立金标逐位一致。
    sets_singleton = (
        "set_id\trole\tgene\n"
        + "".join(f"SSET\tS\t{g}\n" for g in S_SET)
        + "TSET\tT\tT3\n"
    )
    out_tsv, out_log, result = run_node(
        tmp_path, EDGES_TSV, sets_singleton,
        extra_env={"NP_MEASURE": "separation"},
    )
    assert result.returncode == 0, result.stderr  # v1 在此 TypeError 崩溃
    _, rows = read_table(out_tsv)
    assert len(rows) == 1
    row = rows[0]
    assert row["n_s"] == str(len(S_IN)) and row["n_t"] == "1"
    expected = separation_golden(GOLDEN_DIST, S_IN, T_SINGLETON)
    assert expected is not None  # toolbox 0 代入后金标可估
    assert abs(float(row["d_observed"]) - expected) < 1e-12
    assert "not estimable" not in out_log.read_text()


def test_all_na_still_exits_2_with_evidence_written(tmp_path):
    # v3（fix-review R06 失败保行）：集合经最大分量过滤后为空 → NA 行 + log
    # 记因，证据文件照写；但该表没有任何可估组合 → 退出 2（v2 退 2 但无
    # 行级证据；静默 skip 是被指控的缺陷）。
    sets_offgraph = "set_id\trole\tgene\nSSET\tS\tS1\nTSET\tT\tISO1\n"
    out_tsv, out_log, result = run_node(
        tmp_path, EDGES_TSV, sets_offgraph, expect_ok=False,
    )
    assert result.returncode == 2
    _, rows = read_table(out_tsv)  # 证据文件已写出（失败保行）
    assert len(rows) == 1
    row = rows[0]
    assert row["n_s"] == "1" and row["n_t"] == "0"
    for column in ("d_observed", "d_rand_mean", "d_rand_sd", "mc_se", "z",
                   "p_lower", "p_upper", "n_random_used"):
        assert row[column] in ("NA", "nan", ""), f"{column} should be NA, got {row[column]!r}"
    assert "not estimable" in out_log.read_text()
    assert "all rows NA" in result.stderr


def test_partial_na_keeps_row_and_exits_zero(tmp_path):
    # v3 混合情形：一个可估组合 + 一个空集合组合 → 退 0，两行都在表里
    # （空集合组合 NA + log 记因，可估组合正常统计）。
    sets_mixed = (
        "set_id\trole\tgene\n"
        + "".join(f"SSET\tS\t{g}\n" for g in S_SET)
        + "SOFF\tS\tISO1\n"
        + "".join(f"TSET\tT\t{g}\n" for g in T_SET)
    )
    out_tsv, out_log, result = run_node(tmp_path, EDGES_TSV, sets_mixed)
    assert result.returncode == 0, result.stderr
    _, rows = read_table(out_tsv)
    assert len(rows) == 2
    by_set = {row["set_s"]: row for row in rows}
    assert float(by_set["SSET"]["d_observed"]) == closest_golden(GOLDEN_DIST, S_IN, T_IN)
    assert by_set["SOFF"]["n_s"] == "0"
    assert by_set["SOFF"]["d_observed"] in ("NA", "nan", "")
    assert "SOFF vs TSET" in out_log.read_text()


def test_separation_normal_pair_matches_golden(tmp_path):
    out_tsv, out_log, result = run_node(
        tmp_path, EDGES_TSV, SETS_TSV, extra_env={"NP_MEASURE": "separation"},
    )
    assert result.returncode == 0, result.stderr
    _, rows = read_table(out_tsv)
    expected = separation_golden(GOLDEN_DIST, S_IN, T_IN)
    assert abs(float(rows[0]["d_observed"]) - expected) < 1e-12
    b = int(rows[0]["n_random_used"])
    assert 0 < b <= 100
    # 经验 P 界：(b+1)/(B+1) ∈ [1/(B+1), 1]；mc_se = sd/sqrt(B)。
    for column in ("p_lower", "p_upper"):
        p = float(rows[0][column])
        assert 1.0 / (b + 1) - 1e-12 <= p <= 1.0 + 1e-12
    sd = float(rows[0]["d_rand_sd"])
    assert abs(float(rows[0]["mc_se"]) - sd / (b ** 0.5)) < 1e-12


def test_seed_determinism_byte_identical(tmp_path):
    out1, _, _ = run_node(tmp_path / "r1", EDGES_TSV, SETS_TSV)
    out2, _, _ = run_node(tmp_path / "r2", EDGES_TSV, SETS_TSV)
    assert out1.read_bytes() == out2.read_bytes()


def test_default_n_random_is_10000_and_logged(tmp_path):
    # 不传 NP_N_RANDOM → 默认 10000（v2 变更，log 自锚 + 输出列回读）。
    edges = tmp_path / "edges.tsv"
    sets = tmp_path / "sets.tsv"
    edges.write_text(EDGES_TSV, encoding="utf-8")
    sets.write_text(SETS_TSV, encoding="utf-8")
    full_env = dict(os.environ)
    full_env.pop("NP_N_RANDOM", None)
    full_env.update({
        "AUTONOMICS_INPUT0": str(edges), "AUTONOMICS_INPUT1": str(sets),
        "AUTONOMICS_OUTPUT0": str(tmp_path / "p.tsv"),
        "AUTONOMICS_OUTPUT1": str(tmp_path / "p.log"),
    })
    proc = subprocess.run([PYTHON, str(SCRIPT)], env=full_env,
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    log = (tmp_path / "p.log").read_text(encoding="utf-8")
    assert "n_random=10000" in log
    table = (tmp_path / "p.tsv").read_text(encoding="utf-8")
    assert table.splitlines()[1].split("\t")[-1] == "10000"


def test_output_columns_match_manifest_declaration(tmp_path):
    out_tsv, _, _ = run_node(tmp_path, EDGES_TSV, SETS_TSV)
    header, _ = read_table(out_tsv)
    doc = tomllib.loads(MANIFEST.read_text(encoding="utf-8"))["nodes"][0]["doc"]
    # 精确契约：manifest "输出：proximity.tsv ——" 声明的列 == 实际表头。
    declared = doc.split("输出：proximity.tsv ——", 1)[1].split("+ netprox.log", 1)[0]
    declared_columns = [c.strip() for c in declared.strip().split("/")]
    assert declared_columns == header
    # v1 契约缺口（声明 n_in_component 实际输出 n_random_used）不复存在：
    # 声明行里不得再出现未真实输出的列名。
    assert "n_in_component" not in declared_columns
    assert "n_random_used" in declared_columns and "mc_se" in declared_columns
