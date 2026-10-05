# Applying this workflow to public data

Nothing in this folder is committed except this file. Genotype and summary-statistic files are large, and `.gitignore` blocks them.

## Resources (check current versions and access terms before downloading)

| What | Source |
|---|---|
| T2D GWAS summary statistics | DIAGRAM / DIAMANTE consortium (e.g. the multi-ancestry meta-analysis by Mahajan et al., *Nat Genet* 2022), via the DIAGRAM website or the GWAS Catalog |
| Intermediate-trait summary statistics | MAGIC (glycaemic traits, including insulin secretion indices) and GLGC (lipids: HDL-C, triglycerides) |
| Published T2D and process-specific scores | PGS Catalog (pgscatalog.org); search "type 2 diabetes" |
| Variant cluster definitions | Udler et al., *PLoS Med* 2018 (5 clusters); Smith et al., *Nat Med* 2024 (multi-ancestry clusters) |
| Reference / target genotypes | 1000 Genomes Phase 3 in PLINK format (e.g. from the PLINK 2 resources page) |
| Software | PLINK 2 (clumping, scoring); PRSice-2 (C+T); a tutorial in Choi, Mak & O'Reilly, *Nat Protoc* 2020 |

## Outline

```bash
# 1. QC summary statistics: harmonise alleles to the reference, drop ambiguous A/T and C/G SNPs, MAF > 0.01, INFO > 0.8
# 2. Clump against an ancestry-matched reference panel
plink2 --pfile 1kg_eur --clump t2d_sumstats.txt --clump-p1 5e-8 --clump-r2 0.1 --clump-kb 500 --out t2d_clumped
# 3. Score target genotypes
plink2 --pfile target --score t2d_weights.txt 1 2 3 header cols=+scoresums --out t2d_prs
# 4. Partitioned scores: restrict weights to the variants in each published cluster and score separately
# 5. Evaluate in the target cohort: AUC / Nagelkerke R² for T2D; association of each partitioned score with
#    insulin secretion and lipid traits, adjusted for age, sex and genetic principal components
```

**Ancestry.** Scores derived mainly from European GWAS transfer poorly to African-ancestry populations, because LD patterns and allele frequencies differ. Evaluate performance separately by ancestry group and do not assume it transfers.
