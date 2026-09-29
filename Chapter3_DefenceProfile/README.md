# Chapter 3 — Defence profile of *S. suis*

Run from this folder.

```bash
jupyter nbconvert --to notebook --execute --inplace Ssuis_DefenceAnalysis.ipynb
Rscript defence_network.R
```

| Thesis | Code |
|---|---|
| Fig 3.1, Fig S3.1 | `defence_network.R` (panels combined by hand) |
| Fig 3.2–3.5, Table 3.1 | `Ssuis_DefenceAnalysis.ipynb` |
| Fig 3.6 | `hpc/05`–`06`, drawn in iTOL with `supporting_scripts/itol/` |
| Fig 3.7, Fig S3.2–S3.4 | `hpc/07`–`15`, drawn with ThirdKind |
| Tmn-alpha/beta representatives | `supporting_scripts/find_tmn_cluster_medoids.py` |

`hpc/` has the Kelvin2 jobs: genome download, serotyping, Ptolemaea, MMseqs2, trees and GeneRax.
