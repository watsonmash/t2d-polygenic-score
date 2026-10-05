"""End-to-end demonstration on simulated data.

  1. Simulate discovery and target cohorts with two T2D mechanisms.
  2. GWAS in discovery for T2D and two intermediate traits.
  3. Clumping + thresholding -> global T2D PRS; evaluate AUC in the target cohort.
  4. Cluster lead variants by their intermediate-trait association profiles ->
     beta-cell and insulin-resistance partitioned PRS.
  5. Test mechanism specificity: does the beta-cell PRS track insulin secretion
     but not the lipid trait, and vice versa?

Usage: python src/run_demo.py
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

from prs import clump, partition_variants, score
from simulate import gwas, simulate

OUT = Path("outputs"); OUT.mkdir(exist_ok=True)
PAL = {"beta_cell": "#0072B2", "insulin_resistance": "#D55E00", "global": "#5f5f5c"}
plt.rcParams.update({"savefig.dpi": 200, "savefig.bbox": "tight", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": "#d8d8d4",
                     "font.size": 9.5, "axes.titleweight": "bold", "axes.titlelocation": "left",
                     "legend.frameon": False})


def stamp(fig):
    fig.text(0.99, 0.005, "SIMULATED DATA", ha="right", va="bottom", fontsize=8, color="#C0392B")


def main():
    Gd, Pd, Gt, Pt, snps = simulate()
    print(f"discovery n={len(Pd)}, target n={len(Pt)}, SNPs={len(snps)}, "
          f"T2D prevalence={Pd.t2d.mean():.2f}, causal={snps.causal.sum()}")

    ss = {t: gwas(Gd, Pd[t].values.astype(float)) for t in ["t2d", "insulin_secretion", "aip_like"]}
    t2d = ss["t2d"].assign(snp=snps.snp, mechanism=snps.mechanism)

    # --- global C+T PRS across p-value thresholds (thresholds chosen in a held-out half of target)
    tune, test = np.arange(len(Pt)) % 2 == 0, np.arange(len(Pt)) % 2 == 1
    rows = []
    for p1 in [5e-8, 1e-6, 1e-4, 1e-3, 1e-2]:
        idx = clump(t2d, Gd[:5000], p1=p1)
        if len(idx) == 0:
            continue
        s = score(Gt, idx, t2d.beta.values[idx])
        rows.append({"p threshold": p1, "SNPs": len(idx),
                     "AUC (tuning half)": roc_auc_score(Pt.t2d[tune], s[tune]),
                     "AUC (test half)": roc_auc_score(Pt.t2d[test], s[test])})
    ct = pd.DataFrame(rows)
    best = ct.loc[ct["AUC (tuning half)"].idxmax()]
    print("\nClumping + thresholding\n", ct.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print(f"selected p < {best['p threshold']:g}: test-half AUC {best['AUC (test half)']:.3f}")

    # --- partition genome-wide significant lead variants by mechanism
    lead = clump(t2d, Gd[:5000], p1=5e-8)
    Z = np.column_stack([ss["insulin_secretion"].z.values[lead], ss["aip_like"].z.values[lead]])
    # orient each variant to its T2D risk-raising allele
    Z = Z * np.sign(t2d.beta.values[lead])[:, None]
    lab, W, H = partition_variants(Z, k=2)
    # name factors by their trait loadings: [IS+, AIP+, IS-, AIP-]
    names = {}
    for f in range(2):
        names[f] = "beta_cell" if H[f, 2] > H[f, 1] else "insulin_resistance"
    assigned = np.array([names[l] for l in lab])
    truth = snps.mechanism.values[lead]
    tab = pd.crosstab(pd.Series(truth, name="true mechanism"), pd.Series(assigned, name="assigned cluster"))
    print("\nVariant clustering vs simulated truth\n", tab)

    # --- partitioned PRS in target cohort and specificity test
    res = []
    scores = {"global": score(Gt, lead, t2d.beta.values[lead])}
    for m in ["beta_cell", "insulin_resistance"]:
        sel = lead[assigned == m]
        scores[m] = score(Gt, sel, t2d.beta.values[sel])
    for name, s in scores.items():
        zs = (s - s.mean()) / s.std()
        res.append({"Score": name, "SNPs": len(lead) if name == "global" else int((assigned == name).sum()),
                    "T2D AUC": roc_auc_score(Pt.t2d, zs),
                    "r with insulin secretion": stats.pearsonr(zs, Pt.insulin_secretion)[0],
                    "r with lipid trait (AIP-like)": stats.pearsonr(zs, Pt.aip_like)[0]})
    part = pd.DataFrame(res)
    print("\nPartitioned PRS in target cohort\n", part.round(3).to_string(index=False))

    with pd.ExcelWriter(OUT / "prs_results.xlsx") as xw:
        ct.to_excel(xw, sheet_name="clump_threshold", index=False)
        tab.to_excel(xw, sheet_name="cluster_vs_truth")
        part.round(4).to_excel(xw, sheet_name="partitioned_prs", index=False)

    # --- figures
    fig, ax = plt.subplots(figsize=(7, 3))
    ax.scatter(np.arange(len(t2d)), -np.log10(t2d.p), s=4, c=[PAL.get(m, "#bbbbbb") for m in t2d.mechanism])
    ax.axhline(-np.log10(5e-8), color="#5f5f5c", ls=":")
    ax.set_xlabel("SNP position (simulated)"); ax.set_ylabel("−log10 p")
    ax.set_title("Discovery GWAS for T2D (blue: beta-cell loci, orange: insulin-resistance loci)")
    stamp(fig); fig.savefig(OUT / "f1_manhattan.png"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    for m in ["beta_cell", "insulin_resistance"]:
        k = assigned == m
        ax.scatter(Z[k, 0], Z[k, 1], s=24, color=PAL[m], label=f"assigned: {m.replace('_', ' ')}")
    ax.axhline(0, color="#999", lw=0.8); ax.axvline(0, color="#999", lw=0.8)
    ax.set_xlabel("z, insulin secretion (risk allele)"); ax.set_ylabel("z, lipid trait (risk allele)")
    ax.set_title("Lead T2D variants by trait profile"); ax.legend(fontsize=8)
    stamp(fig); fig.savefig(OUT / "f2_variant_clusters.png"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.6, 3.2))
    x = np.arange(len(part)); w = 0.38
    ax.bar(x - w / 2, part["r with insulin secretion"], w, color=PAL["beta_cell"], label="insulin secretion")
    ax.bar(x + w / 2, part["r with lipid trait (AIP-like)"], w, color=PAL["insulin_resistance"], label="lipid trait")
    ax.axhline(0, color="#555", lw=0.8)
    ax.set_xticks(x, [s.replace("_", " ") for s in part.Score]); ax.set_ylabel("Pearson r in target cohort")
    ax.set_title("Mechanism specificity of partitioned scores"); ax.legend(fontsize=8)
    stamp(fig); fig.savefig(OUT / "f3_specificity.png"); plt.close(fig)


if __name__ == "__main__":
    main()
