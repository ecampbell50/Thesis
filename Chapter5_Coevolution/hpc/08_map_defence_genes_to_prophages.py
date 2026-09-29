#!/usr/bin/env python3
# ======================================================================
# Thesis Chapter 5, Methods: prophage annotation (defence gene location)
# Assigns each defence gene to host chromosome or prophage by overlapping PROKKA and geNomad coordinates (writes genome_defence_by_location.csv).
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/Chapter3_Coevolution/map_defence_to_prophage.py
# ======================================================================
# =============================================================================
# Is each defence gene INSIDE a prophage?  (run on the HPC)
#
# Chain:
#   defence_genes_all.fasta  ">{genome}@{locus_tag} TYPE= SUBTYPE= OUTCOME="
#        -> which proteins are defence + their system
#   Prokka GFF  01_prokka/{genome}/*.gff
#        -> locus_tag -> (contig, start, end, strand)
#   geNomad     04_genomad/{genome}_genomad/{genome}_summary/{genome}_virus_summary.tsv
#        -> viral regions:  provirus  CONTIG|provirus_S_E  coordinates "S-E"
#                           whole     CONTIG               coordinates "NA"
#   overlap -> on_prophage per defence gene
#
# Pure stdlib. Writes results into the Chapter3_Coevolution working dir.
# =============================================================================

import os, re, glob, csv
from collections import defaultdict

# ----------------------------------------------------------------------------
PROKKA_DIR  = "/mnt/scratch2/users/40204129/11_DefenseProfileofSsuis/1_SsuisGenomes/ssuis_defence_project/output/01_prokka"
DEF_FASTA   = os.path.join(PROKKA_DIR, "defence_genes_all.fasta")
GENOMAD_DIR = "/mnt/scratch2/users/40204129/Chapter3_Coevolution/04_genomad"
OUT_DIR     = "/mnt/scratch2/users/40204129/Chapter3_Coevolution/05_defence_location"
os.makedirs(OUT_DIR, exist_ok=True)

HDR = re.compile(r"^>(?P<genome>[^@]+)@(?P<locus>\S+)\s+"
                 r"TYPE=(?P<t>\S+)\s+SUBTYPE=(?P<s>\S+)\s+OUTCOME=(?P<o>\S+)")


def parse_defence_fasta(path):
    """genome -> {locus_tag: (type, subtype, outcome)}"""
    by_genome = defaultdict(dict)
    n = 0
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            m = HDR.match(line.strip())
            if not m:                                  # tolerate header variants
                parts = line[1:].split(None, 1)
                gid, locus = parts[0].split("@", 1)
                attrs = dict(re.findall(r"(\w+)=(\S+)", parts[1] if len(parts) > 1 else ""))
                by_genome[gid][locus] = (attrs.get("TYPE", "NA"),
                                         attrs.get("SUBTYPE", "NA"), attrs.get("OUTCOME", "NA"))
            else:
                by_genome[m["genome"]][m["locus"]] = (m["t"], m["s"], m["o"])
            n += 1
    print(f"  defence proteins in fasta: {n}  across {len(by_genome)} genomes")
    return by_genome


def parse_gff(path):
    """locus_tag -> (contig, start, end, strand). Stops at ##FASTA. Returns all contigs too."""
    loci, contigs = {}, set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or "\t" not in line:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            contig = f[0]; contigs.add(contig)
            attrs = f[8]
            m = re.search(r"locus_tag=([^;]+)", attrs) or re.search(r"ID=([^;]+)", attrs)
            if not m:
                continue
            loci[m.group(1)] = (contig, int(f[3]), int(f[4]), f[6])
    return loci, contigs


def parse_genomad(path):
    """Returns (provirus_regions: contig->[(s,e)], whole_viral_contigs: set, all viral contigs)."""
    prov = defaultdict(list); whole = set(); seen = set()
    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            seq = row["seq_name"]; coord = (row.get("coordinates") or "NA").strip()
            contig = seq.split("|", 1)[0]
            seen.add(contig)
            if coord and coord != "NA":
                s, e = coord.split("-")
                prov[contig].append((int(s), int(e)))
            else:
                whole.add(contig)
    return prov, whole, seen


# ----------------------------------------------------------------------------
print("Parsing defence fasta ...")
defence = parse_defence_fasta(DEF_FASTA)

rows = []
n_no_gff = n_no_genomad = 0
n_locus_missing = 0
contig_mismatch_genomes = []          # geNomad viral contig absent from GFF -> naming problem

for genome, loci_sys in sorted(defence.items()):
    gff_hits = glob.glob(os.path.join(PROKKA_DIR, genome, "*.gff*"))
    if not gff_hits:
        n_no_gff += 1
        continue
    loci, gff_contigs = parse_gff(gff_hits[0])

    gpat = glob.glob(os.path.join(GENOMAD_DIR, f"{genome}_genomad", "*_summary",
                                  "*_virus_summary.tsv"))
    if gpat:
        prov, whole, viral_contigs = parse_genomad(gpat[0])
        # contig-name sanity: are geNomad's viral contigs found among GFF contigs?
        if viral_contigs and not (viral_contigs & gff_contigs):
            contig_mismatch_genomes.append(genome)
    else:
        n_no_genomad += 1
        prov, whole = {}, set()

    for locus, (typ, sub, out) in loci_sys.items():
        if locus not in loci:
            n_locus_missing += 1
            continue
        contig, gs, ge, strand = loci[locus]
        on, klass, region = False, "host", ""
        if contig in whole:
            on, klass, region = True, "whole_viral_contig", f"{contig}:whole"
        else:
            for (ps, pe) in prov.get(contig, []):
                if gs <= pe and ge >= ps:                 # overlap
                    on, klass, region = True, "provirus", f"{contig}:{ps}-{pe}"
                    break
        rows.append([genome, locus, typ, sub, out, contig, gs, ge, strand,
                     int(on), klass, region])

# ----------------------------------------------------------------------------
out_genes = os.path.join(OUT_DIR, "defence_gene_prophage_location.csv")
with open(out_genes, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["genome", "locus_tag", "type", "subtype", "outcome", "contig",
                "start", "end", "strand", "on_prophage", "prophage_class", "region"])
    w.writerows(rows)

# per-system summary
agg = defaultdict(lambda: [0, 0])     # (type,subtype) -> [n_total, n_on_prophage]
for r in rows:
    key = (r[2], r[3]); agg[key][0] += 1; agg[key][1] += r[9]
out_sys = os.path.join(OUT_DIR, "defence_system_on_prophage_summary.csv")
with open(out_sys, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["type", "subtype", "n_genes", "n_on_prophage", "frac_on_prophage"])
    for (t, s), (ntot, non) in sorted(agg.items(), key=lambda x: -x[1][1] / max(x[1][0], 1)):
        w.writerow([t, s, ntot, non, round(non / ntot, 4) if ntot else 0])

# ----------------------------------------------------------------------------
tot = len(rows); onp = sum(r[9] for r in rows)
print("\n=== DIAGNOSTICS (check these before trusting the output) ===")
print(f"  genomes processed              : {len(defence) - n_no_gff}")
print(f"  genomes with NO prokka gff     : {n_no_gff}")
print(f"  genomes with NO genomad summary: {n_no_genomad}")
print(f"  defence loci not found in gff  : {n_locus_missing}  (should be ~0)")
print(f"  >>> CONTIG-NAME MISMATCH genomes: {len(contig_mismatch_genomes)} "
      f"(should be 0; geNomad contigs not found in GFF)")
if contig_mismatch_genomes[:5]:
    print("      examples:", contig_mismatch_genomes[:5])
print(f"\n  defence genes placed           : {tot}")
print(f"  on a prophage                  : {onp}  ({100*onp/tot:.1f}%)" if tot else "  no genes")
print("\n  most prophage-associated systems:")
ranked = sorted(agg.items(), key=lambda x: -x[1][1] / max(x[1][0], 1))
for (t, s), (ntot, non) in ranked[:12]:
    if ntot >= 5:
        print(f"    {t}:{s:24s} {non:4d}/{ntot:4d}  ({100*non/ntot:.0f}% on prophage)")
print("\nOutputs ->", OUT_DIR)
print("  - defence_gene_prophage_location.csv")
print("  - defence_system_on_prophage_summary.csv")
