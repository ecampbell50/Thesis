#!/bin/bash
# ======================================================================
# Thesis Chapter 3, Figure 3.7 (marker-gene species tree)
# Concatenates the 40 MAFFT-aligned COGs into one supermatrix with catsequences.
# RECONSTRUCTED: this step was run interactively, not from a script. The command
# below reproduces the outputs found next to it on Kelvin2 (allseqs.fas,
# allseqs.partitions.txt), in .../GeneReconcilliation_of_tmn/Phylogeny_alignments/
# ======================================================================

# 11_catsequences_input_list.txt = the 40 COG*_aligned.faa files (LISTTOCAT.txt)
catsequences 11_catsequences_input_list.txt
