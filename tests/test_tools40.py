import sys
import numpy as np
sys.path.insert(0, "src")
from biomedml.tools40 import ALL_TOOLS, STAT_TOOLS, CLF_TOOLS, EVAL_TOOLS, GRAPH_TOOLS, PANEL_TOOLS

def test_exactly_40():
    assert len(ALL_TOOLS) == 40
    assert len(STAT_TOOLS) == 10 and len(CLF_TOOLS) == 10
    assert len(EVAL_TOOLS) == 10 and len(GRAPH_TOOLS) == 6 and len(PANEL_TOOLS) == 4

def toy():
    rng = np.random.RandomState(0)
    X = rng.randn(60, 20)
    y = np.array([0] * 30 + [1] * 30)
    X[30:, 0] += 1.5  # signal feature
    return X, y

def test_welch_detects_signal():
    X, y = toy()
    t0 = STAT_TOOLS["welch_t"](X[:, 0], y)
    t1 = STAT_TOOLS["welch_t"](X[:, 1], y)
    assert abs(t0) > 3 and abs(t1) < abs(t0)

def test_bh_fdr_basic():
    assert STAT_TOOLS["bh_fdr"]([0.001] * 5 + [0.9] * 95) == 5
    assert STAT_TOOLS["bh_fdr"]([0.5] * 100) == 0

def test_zscore():
    X, _ = toy()
    Z = STAT_TOOLS["zscore"](X)
    assert np.abs(Z.mean(0)).max() < 1e-9

def test_majority_baseline():
    X, y = toy()
    assert CLF_TOOLS["majority"](X, y)[0] == 0.5

def test_eval_perfect_score():
    y = np.array([0, 0, 1, 1])
    s = np.array([0.1, 0.2, 0.8, 0.9])
    assert EVAL_TOOLS["roc_auc"](y, s) == 1.0
    assert EVAL_TOOLS["confusion"](y, s) == [[2, 0], [0, 2]]

def test_graph_tools():
    C = np.eye(10); C[0, 1] = C[1, 0] = 0.9
    assert GRAPH_TOOLS["edge_density"](C, 0.5) > 0
    r = GRAPH_TOOLS["community"](C, 0.5)
    assert r["n_components"] == 9  # one pair + 8 isolates

def test_sign_stability():
    coefs = np.array([[1, -1, 1], [1, -1, -1], [1, -1, 1]], dtype=float)
    stab = PANEL_TOOLS["sign_stability"](coefs)
    assert stab[0] == 1.0 and stab[1] == 1.0
