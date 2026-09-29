# Chapter 2 — Ptolemaea

The pipeline itself: [github.com/ecampbell50/Ptolemaea](https://github.com/ecampbell50/Ptolemaea)

- `hpc/` — genome retrieval and the per-species Ptolemaea runs
- `scripts/` — statistics and Figure 2.1 panels, from the per-species outputs in `data/`

```bash
python3 scripts/analyse_results.py            # stats -> results_stats.txt, figures/
python3 scripts/build_network.py --figure     # network panel -> figures/
```
