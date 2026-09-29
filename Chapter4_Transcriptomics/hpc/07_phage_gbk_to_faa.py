#!/usr/bin/env python3
# ======================================================================
# Thesis Chapter 4, Figures S4.1, S4.2
# Converts the Pharokka GenBank output to protein FASTA (for eggNOG).
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/3_Transcriptomics/7_GenomeAnnotations/2_PhageGenomes/1_Pharokka/gbk_to_faa.py
# ======================================================================
from Bio import SeqIO
import sys

if len(sys.argv) != 3:
    print("Usage: python gbk_to_faa.py input.gbk output.faa")
    sys.exit(1)

gbk_file = sys.argv[1]
faa_file = sys.argv[2]

count = 0
with open(faa_file, 'w') as output:
    for record in SeqIO.parse(gbk_file, 'genbank'):
        for feature in record.features:
            if feature.type == 'CDS':
                if 'translation' in feature.qualifiers:
                    locus_tag = feature.qualifiers.get('locus_tag', ['Unknown'])[0]
                    product = feature.qualifiers.get('product', ['hypothetical protein'])[0]
                    translation = feature.qualifiers['translation'][0]
                    output.write(f">{locus_tag} {product}\n{translation}\n")
                    count += 1

print(f"Extracted {count} proteins to {faa_file}")
