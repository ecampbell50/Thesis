# Stats cheat-sheet — *S. suis* prophage/defence project

A plain-English reference for the tests used in this analysis. Two questions to ask of *every*
result: **(1) Is it real?** → the p-value (or q-value). **(2) Is it big enough to matter?** →
the effect size. With ~2,100 genomes, even tiny effects get tiny p-values, so always read both.

---

## The two foundations

| Concept | What it is | Rule of thumb |
|---|---|---|
| **p-value** | Probability of seeing a pattern this strong *by pure luck* if nothing were really going on. | < 0.05 = "luck is an unconvincing explanation." `2e-44` = absurdly unlikely to be luck. **Real ≠ big.** |
| **q-value (BH-FDR)** | A p-value *corrected* for the fact that you ran many tests at once (testing 100 things, ~5 look "significant" by chance). | Use the q-value whenever many things were tested (e.g. 100+ defence systems). q < 0.05 = survives the correction. |

---

## Family 1 — "Are these GROUPS different?"

| Test | Use it when… | How it works (1 line) | Effect size | We used it for |
|---|---|---|---|---|
| **Mann-Whitney U** | comparing **2 groups** of numbers | ranks all values; do one group's sit higher? | difference in medians | 1 vs >1 prophage: defence count |
| **Kruskal-Wallis** | comparing **3+ groups** of numbers | same as M-W but for many groups; "is *anyone* different?" | **η² (eta-squared)** = % of variation explained by the grouping | 0/1/>1 burden; 14 serotypes; prophage size across communities (η²=0.86!) |
| **Fisher's exact** | **yes/no** trait, **2 groups** | exact odds of a 2×2 split this lopsided | **odds ratio (OR)**: OR=472 → 472× more likely | "is system X a serotype-16 marker?" |
| **Chi-square** | **yes/no** trait, **2+ groups**, large counts | approximation of Fisher for bigger tables | (compare fractions) | system present/absent across 0/1/>1 burden |

> Kruskal-Wallis tells you *that* some group differs, not *which* — follow up with Mann-Whitney on the specific pairs.
> Chi-square breaks with tiny counts (empty cells) — use Fisher there instead.

---

## Family 2 — "Do two NUMBERS move together?"

| Test | Use it when… | How it works | Reading the number |
|---|---|---|---|
| **Spearman ρ ("rho")** | two numeric variables, "do they rise/fall together?" | ranks both, measures how well rankings line up | **−1** opposite · **0** unrelated · **+1** together. Doesn't need a straight line, just a consistent trend. |
| **Partial Spearman** | same, but **removing the effect of a confounder** | strips out the part each variable shares with the confounder, then correlates the leftovers | "holding [confounder] constant, do they *still* move together?" |

> **This pair is the heart of the confounder check.** Burden vs defence = +0.27 (plain Spearman),
> but burden vs defence *holding genome size constant* = −0.20 (partial Spearman). The positive
> link lived entirely in genome size.

---

## Family 3 — "Untangle EVERYTHING at once"

| Test | Use it when… | What it gives you |
|---|---|---|
| **Negative-binomial regression** | outcome is a **count** (defence genes), and you want to weigh **many factors simultaneously** (burden + genome size + serotype + assembly quality…) | the effect of each factor *with the others held fixed* |

**Reading it:**
- **IRR (Incidence Rate Ratio)** = the multiplier on the outcome per 1-unit increase in a predictor, holding everything else fixed.
  - IRR **= 1.00** → no effect · **> 1** → more · **< 1** → less.
  - e.g. burden IRR 0.981 = "each extra prophage → ~2% *fewer* defence genes, at the same genome size."
- **95% confidence interval** = plausible range for that multiplier. **Does it include 1.00?**
  - excludes 1.00 → effect is solid · includes 1.00 → could just be noise.
- *Why "negative-binomial"?* Counts aren't bell-curve data, and defence counts are extra-spread-out ("over-dispersed"). NB is the standard, well-behaved model for messy biological counts. (Plain linear regression assumes smooth measurements — wrong tool for counts.)

---

## Two traps that bit us (worth remembering)

| Trap | What it means | What we did |
|---|---|---|
| **Confounding** | A third variable causes *both* things you're comparing, faking a link. Genome size → more prophages AND more defence, so they correlate without one causing the other. | Control for it (partial Spearman / regression). The link vanished → genome-size artifact. |
| **Pseudoreplication** | Treating non-independent data points as independent fakes having more evidence than you have. Prophages in one community are near-copies (η²=0.86 for size). | Re-test at the **community level** (one point per community), not per-prophage. |

---

## The session in one sentence
Spearman said burden & defence move together → partial Spearman + negative-binomial regression
said *"only because both track genome size"* → confounder-free answer = **no link**.

## Where each was used (scripts in repo root)
- `analyse_burden_defence.py` — Kruskal-Wallis, Mann-Whitney, Spearman, chi-square (+BH-FDR)
- `analyse_serotype_defence.py` — Kruskal-Wallis, Fisher's exact (+BH-FDR)
- `analyse_defence_location.py` — Spearman, chi-square
- `analyse_confounders.py` — partial Spearman, **negative-binomial regression**, Kruskal-Wallis/η², Mann-Whitney
