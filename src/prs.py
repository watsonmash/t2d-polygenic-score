"""Clumping + thresholding PRS, and mechanism-partitioned PRS via clustering of
variant association profiles across intermediate traits (simplified version of the
soft-clustering approach used to define T2D 'process-specific' polygenic scores).
"""
import numpy as np
import pandas as pd
from sklearn.decomposition import NMF


def ld_r2(G, i, j):
    return np.corrcoef(G[:, i], G[:, j])[0, 1] ** 2


def clump(sumstats, G_ref, p1=5e-8, r2=0.1, window=50):
    """Greedy LD clumping: take the most significant SNP, drop neighbours in LD."""
    order = np.argsort(sumstats["p"].values)
    removed = np.zeros(len(sumstats), bool)
    index = []
    for i in order:
        if sumstats["p"].iat[i] > p1:
            break
        if removed[i]:
            continue
        index.append(i)
        lo, hi = max(0, i - window), min(len(sumstats), i + window + 1)
        for j in range(lo, hi):
            if j != i and not removed[j] and ld_r2(G_ref, i, j) > r2:
                removed[j] = True
    return np.array(index)


def score(G, idx, weights):
    return G[:, idx].astype(float) @ weights


def partition_variants(z_matrix, k=2, seed=2026):
    """Non-negative factorisation of variant x trait z-scores (each trait split into
    positive and negative channels, as in bNMF). Variants are assigned to the factor
    with the largest weight; trait loadings label the factors."""
    Z = np.asarray(z_matrix, float)
    Zp = np.hstack([np.clip(Z, 0, None), np.clip(-Z, 0, None)])
    model = NMF(n_components=k, init="nndsvda", random_state=seed, max_iter=2000)
    W = model.fit_transform(Zp)
    H = model.components_
    return W.argmax(1), W, H
