"""Differential expression and multiple-testing control."""
from __future__ import annotations
import numpy as np
from scipy import stats as sps


def de_genes(X, y, max_genes=None):
    """Welch t-test per gene; returns arrays (t, p) aligned to genes."""
    X0, X1 = X[y == 0], X[y == 1]
    t, p = sps.ttest_ind(X1, X0, axis=0, equal_var=False, nan_policy="omit")
    return np.nan_to_num(t), np.nan_to_num(p, nan=1.0)


def bh_fdr(p, alpha=0.05):
    """Benjamini-Hochberg; returns boolean mask of discoveries."""
    p = np.asarray(p)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    thresh = alpha * np.arange(1, n + 1) / n
    below = ranked <= thresh
    if not below.any():
        return np.zeros(n, bool)
    kmax = np.max(np.where(below))
    mask = np.zeros(n, bool)
    mask[order[:kmax + 1]] = True
    return mask


def coexpression_graph(X, top_idx, thresh=0.6):
    """Gene co-expression adjacency over selected genes (|r| >= thresh)."""
    Z = X[:, top_idx]
    Z = (Z - Z.mean(0)) / (Z.std(0) + 1e-8)
    R = (Z.T @ Z) / Z.shape[0]
    A = (np.abs(R) >= thresh).astype(np.float32)
    np.fill_diagonal(A, 0.0)
    return A, R
