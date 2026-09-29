# Thesis Code Repository

Analysis code for the PhD thesis:

> **Anti-phage defence, host infection response, and prophage co-evolutionary dynamics in *Streptococcus suis***
> Emmet B. T. Campbell, Queen's University Belfast, 2026

| Folder | Thesis chapter |
|---|---|
| `Chapter2_Ptolemaea/` | Ptolemaea on 700 ESKAPE + *E. coli* genomes (pipeline: [github.com/ecampbell50/Ptolemaea](https://github.com/ecampbell50/Ptolemaea)) |
| `Chapter3_DefenceProfile/` | Defence profile of 2,119 *S. suis* genomes |
| `Chapter4_Transcriptomics/` | Host and phage transcriptomics |
| `Chapter5_Coevolution/` | Prophage–host co-evolution |

In each chapter, `hpc/` holds the SLURM jobs run on the Kelvin2 HPC, numbered in the
order they were run. Every script starts with a short header naming the thesis figure
or table it feeds. Generated figures and tables are not committed; re-run the scripts
to produce them.

## License

See [LICENSE](LICENSE).
