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
    d_ab = closest_golden(dist, s_members, t_members)
    d_aa = closest_golden(dist, s_members, s_members, exclude_self=True)
    d_bb = closest_golden(dist, t_members, t_members, exclude_self=True)
    if d_ab is None or d_aa is None or d_bb is None:
        return None
    return d_ab - 0.5 * (d_aa + d_bb)


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


def test_separation_matches_golden_and_singleton_is_none():
    mod = load_module()
    import networkx as nx

    graph = nx.Graph()
    graph.add_edges_from(GOLDEN_EDGES)
    component = max(nx.connected_components(graph), key=len)
    graph = graph.subgraph(component).copy()
    expected = separation_golden(GOLDEN_DIST, S_IN, T_IN)
    value, reason = mod.proximity(graph, "separation", S_IN, T_IN)
    assert value == expected and reason is None
    # 单元素集合：排除自身后集合内距离无定义 → (None, "singleton")，不崩溃。
    value, reason = mod.proximity(graph, "separation", S_IN, T_SINGLETON)
    assert value is None and reason == "singleton"


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


def test_separation_singleton_outputs_na_row_without_crash(tmp_path):
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
    for column in ("d_observed", "d_rand_mean", "d_rand_sd", "mc_se", "z",
                   "p_lower", "p_upper", "n_random_used"):
        assert row[column] in ("NA", "nan", ""), f"{column} should be NA, got {row[column]!r}"
    assert "not estimable" in out_log.read_text()


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


def test_empty_after_component_filter_exits(tmp_path):
    sets_offgraph = "set_id\trole\tgene\nSSET\tS\tS1\nTSET\tT\tISO1\n"
    _, _, result = run_node(
        tmp_path, EDGES_TSV, sets_offgraph, expect_ok=False,
    )
    # ISO1 在小组件里被分量过滤 → 该组合 skipped；无组合可估 → 退出 2。
    assert result.returncode == 2
