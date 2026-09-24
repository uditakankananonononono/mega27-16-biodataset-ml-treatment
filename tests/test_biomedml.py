import numpy as np
from biomedml.stats import de_genes, bh_fdr, coexpression_graph
from biomedml.models import ExprCNN, LogRegNP, train_torch, predict_torch, auc_score


def toy(seed=0, n=120, g=200, n_signal=12, latent_scale=0.0):
    rng = np.random.default_rng(seed)
    y = np.repeat([0, 1], n // 2)
    X = rng.normal(0, 1, (n, g)).astype(np.float32)
    X[y == 1, :n_signal] += 1.5  # planted DE genes
    if latent_scale:
        latent = rng.normal(0, 1, n)  # shared factor -> correlated block
        X[:, :n_signal] += latent_scale * latent[:, None]
    return X, y


def test_de_genes_finds_planted():
    X, y = toy()
    t, p = de_genes(X, y)
    top = set(np.argsort(p)[:20])
    assert len(top & set(range(12))) >= 8


def test_bh_fdr_controls_and_discovers():
    X, y = toy()
    _, p = de_genes(X, y)
    mask = bh_fdr(p, 0.05)
    assert mask[:12].all()
    null_p = np.random.default_rng(1).uniform(0, 1, 1000)
    assert bh_fdr(null_p, 0.05).sum() <= 20  # generous FPR bound


def test_coexpression_graph_recovers_block():
    X, y = toy(latent_scale=0.9)
    A, R = coexpression_graph(X, np.arange(12), thresh=0.5)
    assert A.sum() > 20  # planted correlated block connects


def test_cnn_and_logreg_learn():
    X, y = toy(n=160)
    tr = np.arange(120); te = np.arange(120, 160)
    cnn = ExprCNN(X.shape[1])
    train_torch(cnn, X, y, tr, epochs=40, lr=3e-3)
    acc_cnn = (predict_torch(cnn, X[te]) == y[te]).mean()
    lr = LogRegNP(X.shape[1], epochs=400)
    lr.fit(X[tr], y[tr])
    acc_lr = (lr.predict(X[te]) == y[te]).mean()
    assert acc_cnn > 0.85 and acc_lr > 0.85


def test_auc_perfect_and_chance():
    assert abs(auc_score([0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9]) - 1.0) < 1e-9
    a = auc_score([0, 1] * 50, [0.5] * 100)
    assert abs(a - 0.5) < 1e-9
