# Chapter 2 — Ptolemaea

The Ptolemaea defence-system annotation pipeline described in Chapter 2 is maintained
as its own repository, with documentation, a Singularity container and SLURM submission
scripts:

**https://github.com/ecampbell50/Ptolemaea**

It is kept separate because it is a general-purpose tool intended for reuse beyond this
thesis. Chapter 3 of this repository contains the *output* of running it across the
2,119-genome *S. suis* collection
(`Chapter3_DefenceProfile/data/Ssuis_rawptolemaeaoutput.csv`, 97,099 annotated proteins).

## Upstream SLURM pipeline

The numbered SLURM scripts that produced the Chapter 3 annotations (genome download,
PROKKA, BBMap, serotyping, PADLOC, DefenseFinder, bidirectional BLAST, Roary, FastTree,
Coinfinder) currently live in `github.com/ecampbell50/Ptolemaea_old` under
`scripts/slurm/` and `scripts/bash/`. They are not yet in either current repository.
