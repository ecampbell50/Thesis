# Supporting scripts

Methods-supporting code that is *not* a main figure. Extracted from the original
working notebook and kept here for reproducibility.

| Script | What it does | Methods link |
|--------|--------------|--------------|
| `find_tmn_cluster_medoids.py` | Finds the medoid (most representative) sequence of the large and small CLANS *tmn* clusters via pairwise identity. Returns the *tmn*-alpha / *tmn*-beta representatives (`1005041.3@…00375`, `1214195.3@…00646`). | "Two structural variants of *Tmn*" / GeneRax reconciliation |
| `tmn_defence_comparison.py` | Exploratory comparison of defence content in *tmn*-alpha vs *tmn*-beta genomes. **Superseded** in the manuscript by the within-genome-size-bin Mann–Whitney test in the notebook (§5); kept for reference. | Fig 3 section |
| `itol/make_itol_external_shapes.py` | Generates an iToL `DATASET_EXTERNALSHAPE` file (presence/absence of selected subtypes) to decorate the Figure 6 tree. | Fig 6 |
| `itol/make_itol_serotype_colorstrip.py` | Generates an iToL `DATASET_COLORSTRIP` file colouring tree leaves by serotype. | Fig 6 |
| `itol/annotations/` | The actual iToL annotation files used on the Figure 6 tree (serotype, genome size, *tmn* variant, pathogenicity). | Fig 6 |

These depend on the same dataframes built in the main notebook (run the notebook's
load section first, or adapt the paths). The *medoid* script additionally needs
Biopython and the per-cluster FASTA / ID lists.

## Figure 6 (phylogenetic tree) and Figure 7 (reconciliation)

These figures are produced by external tools, not by code here:

- **Figure 6** — FastTree core-gene phylogeny visualised in **iToL**, decorated with
  the annotation files in `itol/annotations/`. The tree file itself is built by
  Roary + FastTree (see Methods, "S. suis core-genome phylogeny").
- **Figure 7** — GeneRax gene-tree/species-tree reconciliation, rendered with
  **ThirdKind**, with AlphaFold structures inset (see Methods, "GeneRax … reconciliation").

### iToL annotation files

The `annotations/` folder contains several versions accumulated during drafting
(e.g. `FINAL_TmnGroups_All3.txt`, `FINALFINAL_TmnGroups_All3.txt`,
`EDITED_TmnGroups_iToLDataset.txt`). The `FINALFINAL_*`/`EDITED_*` files are the
latest; the others are kept so nothing is lost. Confirm which you uploaded to
iToL before relying on one.
