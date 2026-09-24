"""1D-CNN over expression profiles + logistic baseline (numpy)."""
from __future__ import annotations
import numpy as np
import torch
import torch.nn as nn


class ExprCNN(nn.Module):
    def __init__(self, n_genes, n_classes=2, hidden=32):
        super().__init__()
        # treat the (sorted) expression vector as a 1-channel signal
        self.conv = nn.Sequential(
            nn.Conv1d(1, hidden, 25, stride=5, padding=12), nn.ReLU(),
            nn.Conv1d(hidden, hidden, 15, stride=5, padding=7), nn.ReLU(),
            nn.AdaptiveAvgPool1d(8), nn.Flatten())
        self.head = nn.Linear(hidden * 8, n_classes)

    def forward(self, X):  # (B, G)
        return self.head(self.conv(X.unsqueeze(1)))


def train_torch(model, X, y, tr, epochs=20, lr=1e-3, batch=32, seed=0):
    torch.manual_seed(seed)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    lossf = nn.CrossEntropyLoss()
    Xt = torch.tensor(X, dtype=torch.float32)
    yt = torch.tensor(y, dtype=torch.long)
    for _ in range(epochs):
        perm = torch.randperm(len(tr))
        for i in range(0, len(tr), batch):
            b = torch.tensor(tr)[perm[i:i + batch]]
            opt.zero_grad()
            loss = lossf(model(Xt[b]), yt[b])
            loss.backward()
            opt.step()
    return model


def predict_torch(model, X):
    model.eval()
    with torch.no_grad():
        Xt = torch.tensor(X, dtype=torch.float32)
        return model(Xt).argmax(1).numpy()


class LogRegNP:
    def __init__(self, n_in, lr=0.1, epochs=500, l2=1e-3, seed=0):
        self.W = np.random.default_rng(seed).normal(0, 0.01, (n_in, 2))
        self.b = np.zeros(2)
        self.lr, self.epochs, self.l2 = lr, epochs, l2

    def fit(self, X, y):
        Y = np.eye(2)[y]
        for _ in range(self.epochs):
            Z = X @ self.W + self.b
            Z -= Z.max(1, keepdims=True)
            P = np.exp(Z) / np.exp(Z).sum(1, keepdims=True)
            G = (P - Y) / len(X)
            self.W -= self.lr * (X.T @ G + self.l2 * self.W)
            self.b -= self.lr * G.sum(0)
        return self

    def predict(self, X):
        return (X @ self.W + self.b).argmax(1)


def auc_score(y_true, scores):
    """Rank-based AUC (Mann-Whitney)."""
    y_true = np.asarray(y_true); scores = np.asarray(scores)
    pos = scores[y_true == 1]; neg = scores[y_true == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    wins = sum((p > neg).sum() + 0.5 * (p == neg).sum() for p in pos)
    return float(wins / (len(pos) * len(neg)))
