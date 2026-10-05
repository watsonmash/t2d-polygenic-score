"""Simulate genotypes, two biological mechanisms of type 2 diabetes, and GWAS summary statistics.

Design (deliberately simple and fully transparent):
  * M SNPs in LD blocks; genotypes from thresholded correlated Gaussians (0/1/2).
  * Causal SNPs act through one of two mechanisms:
      - beta-cell dysfunction  -> lowers an insulin-secretion trait
      - insulin resistance     -> raises an adiposity/lipid trait (AIP-like)
  * T2D liability = both mechanism traits + noise; prevalence set by threshold.
  * A discovery cohort supplies GWAS summary statistics for T2D and for both
    intermediate traits; an independent target cohort is used for scoring.
"""
import numpy as np
import pandas as pd
from scipy.stats import norm


def simulate_genotypes(n, m, rng, block=20, rho=0.6, maf=None):
    maf = rng.uniform(0.05, 0.5, m) if maf is None else maf
    G = np.empty((n, m), dtype=np.int8)
    for s in range(0, m, block):
        k = min(block, m - s)
        idx = np.arange(k)
        C = rho ** np.abs(idx[:, None] - idx[None, :])          # AR(1) LD within block
        L = np.linalg.cholesky(C)
        for hap in range(2):
            Z = rng.standard_normal((n, k)) @ L.T
            H = (Z < norm.ppf(maf[s:s + k])).astype(np.int8)
            G[:, s:s + k] = H if hap == 0 else G[:, s:s + k] + H
    return G, maf


def simulate(n_disc=20000, n_target=5000, m=4000, n_causal=(40, 40), h2=(0.25, 0.25),
             prevalence=0.10, seed=2026):
    rng = np.random.default_rng(seed)
    G_all, maf = simulate_genotypes(n_disc + n_target, m, rng)
    blocks = np.arange(m) // 20
    # choose causal SNPs in distinct blocks so each mechanism has its own loci
    cb = rng.choice(np.unique(blocks), sum(n_causal), replace=False)
    causal = np.array([rng.choice(np.where(blocks == b)[0]) for b in cb])
    mech = np.array(["beta_cell"] * n_causal[0] + ["insulin_resistance"] * n_causal[1])

    X = (G_all - 2 * maf) / np.sqrt(2 * maf * (1 - maf))           # standardised genotypes
    traits = {}
    for j, (name, sign) in enumerate([("beta_cell", -1), ("insulin_resistance", +1)]):
        sel = causal[mech == name]
        beta = rng.normal(0, np.sqrt(h2[j] / len(sel)), len(sel))
        g = X[:, sel] @ beta
        e = rng.normal(0, np.sqrt(1 - g.var()), len(g))
        traits[name] = g + e
        traits[f"g_{name}"] = g
    insulin_secretion = -traits["beta_cell"]                       # low secretion = dysfunction
    aip_like = traits["insulin_resistance"]
    liability = 0.7 * traits["beta_cell"] + 0.7 * traits["insulin_resistance"] + rng.normal(0, 0.6, len(aip_like))
    t2d = (liability > np.quantile(liability, 1 - prevalence)).astype(int)

    pheno = pd.DataFrame({"t2d": t2d, "insulin_secretion": insulin_secretion, "aip_like": aip_like,
                          "true_g_beta_cell": traits["g_beta_cell"],
                          "true_g_insulin_resistance": traits["g_insulin_resistance"]})
    snps = pd.DataFrame({"snp": [f"rs{100000 + i}" for i in range(m)], "block": blocks, "maf": maf,
                         "causal": False, "mechanism": "none"})
    snps.loc[causal, "causal"] = True
    snps.loc[causal, "mechanism"] = mech
    disc, targ = slice(0, n_disc), slice(n_disc, None)
    return (G_all[disc], pheno.iloc[disc].reset_index(drop=True),
            G_all[targ], pheno.iloc[targ].reset_index(drop=True), snps)


def gwas(G, y, binary=False):
    """Marginal per-SNP association: linear regression on standardised genotype.
    For a binary outcome this is a linear probability model; the z-statistics are
    what clumping and clustering use, and they rank SNPs as logistic GWAS would."""
    Xs = (G - G.mean(0)) / G.std(0)
    yc = y - y.mean()
    n = len(y)
    beta = Xs.T @ yc / n
    resid_var = ((yc[:, None] - Xs * beta) ** 2).sum(0) / (n - 2)
    se = np.sqrt(resid_var / n)
    z = beta / se
    p = 2 * norm.sf(np.abs(z))
    return pd.DataFrame({"beta": beta, "se": se, "z": z, "p": p})
