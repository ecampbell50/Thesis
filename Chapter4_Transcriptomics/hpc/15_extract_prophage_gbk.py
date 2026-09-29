#!/usr/bin/env python3
# ======================================================================
# Thesis Chapter 4, Figure S4.3
# Cuts prophage loci out of the strain GenBank files for clinker.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/3_Transcriptomics/7_GenomeAnnotations/extract_prophage_gbk.py
# ======================================================================
# Run from: /mnt/scratch2/users/40204129/3_Transcriptomics/7_GenomeAnnotations/

from Bio import SeqIO

BASE = "/mnt/scratch2/users/40204129/3_Transcriptomics/7_GenomeAnnotations"

prophage_ids = {
    "C67": {"C67_prot_00947","C67_prot_00948","C67_prot_00949",
            "C67_prot_00950","C67_prot_00951","C67_prot_00952","C67_prot_00953",
            "C67_prot_00954","C67_prot_00955","C67_prot_00956","C67_prot_00957",
            "C67_prot_00958","C67_prot_00959","C67_prot_00960","C67_prot_00961",
            "C67_prot_00962","C67_prot_00963","C67_prot_00964","C67_prot_00965"},
    "D32": {"D32_prot_01483","D32_prot_01484","D32_prot_01485","D32_prot_01486",
            "D32_prot_01487","D32_prot_01488","D32_prot_01489","D32_prot_01490",
            "D32_prot_01491","D32_prot_01492","D32_prot_01493","D32_prot_01494",
            "D32_prot_01495","D32_prot_01496","D32_prot_01497","D32_prot_01498",
            "D32_prot_01499","D32_prot_01500","D32_prot_01501"},
    "D68": {"D68_prot_01332","D68_prot_01333","D68_prot_01334","D68_prot_01335",
            "D68_prot_01336","D68_prot_01337","D68_prot_01338","D68_prot_01339",
            "D68_prot_01340","D68_prot_01341","D68_prot_01342","D68_prot_01343",
            "D68_prot_01344","D68_prot_01345","D68_prot_01346","D68_prot_01347",
            "D68_prot_01348","D68_prot_01349","D68_prot_01350"},
}

gbff = {"C67": f"{BASE}/1_StrainGenomes/C67_bakta/C67.gbff",
        "D32": f"{BASE}/1_StrainGenomes/D32_bakta/D32.gbff",
        "D68": f"{BASE}/1_StrainGenomes/D68_bakta/D68.gbff"}

for strain, ids in prophage_ids.items():
    for rec in SeqIO.parse(gbff[strain], "genbank"):
        loci = [f for f in rec.features
                if f.type == "CDS" and f.qualifiers.get("locus_tag", [""])[0] in ids]
        if not loci:
            continue
        start = min(int(f.location.start) for f in loci)
        end   = max(int(f.location.end)   for f in loci)
        sliced = rec[start:end]
        sliced.id = f"{strain}_prophage"
        sliced.name = f"{strain}_prophage"
        out = f"{BASE}/3_ProphageComparison/{strain}_prophage.gbk"
        SeqIO.write(sliced, out, "genbank")
        print(f"{strain}: {rec.id}:{start}-{end} → {out}")
