#!/usr/bin/env python3
# ======================================================================
# Thesis Chapter 4, Table S4.3
# Assigns Bakta genes to geNomad prophage regions.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/3_Transcriptomics/7_GenomeAnnotations/3_ConsolidatingAnnotations/C67_consolidation/match_genomad_to_bakta.py
# ======================================================================
"""
Cross-reference geNomad viral gene coordinates with Bakta GFF3 annotations.

geNomad gene coordinates for provirus genes are absolute positions within
their parent contig. This script matches them to Bakta CDS features by
coordinate overlap.

Usage:
    python match_genomad_to_bakta.py \
        --gff3 C67_annotated.gff3 \
        --genomad C67_HBHC_prefixed_virus_genes.tsv \
        --output C67_viral_genes_annotated.tsv
"""

import argparse
import csv
import sys
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class BaktaFeature:
    contig: str
    start: int       # 1-based, inclusive (GFF3 coords)
    end: int         # 1-based, inclusive
    strand: str
    locus_tag: str
    product: str
    raw_attrs: str


@dataclass
class GenomadGene:
    gene_id: str
    parent_contig: str       # e.g. "contig_4.length_152347"
    provirus_region: Optional[str]  # e.g. "provirus_55949_69028" or None
    start: int               # absolute coords within parent contig (1-based)
    end: int
    strand: int              # 1 or -1
    marker: str
    evalue: str
    bitscore: str
    virus_hallmark: str
    annotation_description: str


def parse_gff3(path: str) -> tuple[dict[str, list[BaktaFeature]], dict[str, int]]:
    """
    Parse GFF3 and return:
      - dict of CDS features keyed by contig name
      - dict of contig name -> length from ##sequence-region headers
    """
    features: dict[str, list[BaktaFeature]] = {}
    contig_lengths: dict[str, int] = {}

    with open(path) as f:
        for line in f:
            # Capture ##sequence-region for exact contig lengths
            if line.startswith("##sequence-region"):
                parts = line.strip().split()
                if len(parts) == 4:
                    contig_lengths[parts[1]] = int(parts[3])
                continue
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "CDS":
                continue

            contig = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attrs = parts[8]

            # Extract locus_tag and product from attributes
            locus_tag = ""
            product = ""
            for attr in attrs.split(";"):
                if attr.startswith("locus_tag="):
                    locus_tag = attr[len("locus_tag="):]
                elif attr.startswith("product="):
                    product = attr[len("product="):]

            feat = BaktaFeature(
                contig=contig,
                start=start,
                end=end,
                strand=strand,
                locus_tag=locus_tag,
                product=product,
                raw_attrs=attrs,
            )
            features.setdefault(contig, []).append(feat)

    return features, contig_lengths


def parse_genomad(path: str) -> list[GenomadGene]:
    """Parse geNomad virus_genes.tsv."""
    genes = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene_id = row["gene"]

            # Parse parent contig and optional provirus region from gene ID
            # Format examples:
            #   C67_HBHC|21171_DNR38.contig_4.length_152347|provirus_55949_69028_64
            #   C67_HBHC|21171_DNR38.contig_46.length_431_1
            #
            # Strategy: split on '|', take parts after the sample prefix
            parts = gene_id.split("|")
            # parts[0] = sample prefix e.g. "C67_HBHC"
            # parts[1] = "21171_DNR38.contig_4.length_152347"  OR
            #            "21171_DNR38.contig_46.length_431"
            # parts[2] (if present) = "provirus_55949_69028_64"

            contig_part = parts[1] if len(parts) > 1 else gene_id
            provirus_region = None
            gene_suffix = None

            if len(parts) == 3:
                # Has provirus region
                prov_and_num = parts[2]  # e.g. "provirus_55949_69028_64"
                # Split off gene number from end
                prov_tokens = prov_and_num.rsplit("_", 1)
                provirus_region = prov_tokens[0]  # "provirus_55949_69028"
                gene_suffix = prov_tokens[1]      # "64"
            else:
                # No provirus - gene number is suffix of contig_part
                # e.g. "21171_DNR38.contig_46.length_431_1" -> gene num = "1"
                pass

            # Extract just the contig name for matching against GFF3
            # contig_part looks like "21171_DNR38.contig_4.length_152347"
            # We want "contig_4" style name - but Bakta may name it differently.
            # We'll store the full contig_part so the user can match manually if needed.
            # Also extract a short contig label.
            contig_label = contig_part  # keep full for now

            genes.append(GenomadGene(
                gene_id=gene_id,
                parent_contig=contig_label,
                provirus_region=provirus_region,
                start=int(row["start"]),
                end=int(row["end"]),
                strand=int(row["strand"]),
                marker=row.get("marker", "NA"),
                evalue=row.get("evalue", "NA"),
                bitscore=row.get("bitscore", "NA"),
                virus_hallmark=row.get("virus_hallmark", "NA"),
                annotation_description=row.get("annotation_description", "NA"),
            ))
    return genes


def overlaps(g_start: int, g_end: int, b_start: int, b_end: int,
             min_overlap_frac: float = 0.5) -> bool:
    """
    Check if geNomad gene [g_start, g_end] overlaps Bakta feature [b_start, b_end]
    by at least min_overlap_frac of the shorter feature.
    Coordinates are 1-based inclusive.
    """
    overlap_start = max(g_start, b_start)
    overlap_end = min(g_end, b_end)
    if overlap_end < overlap_start:
        return False
    overlap_len = overlap_end - overlap_start + 1
    shorter = min(g_end - g_start + 1, b_end - b_start + 1)
    return (overlap_len / shorter) >= min_overlap_frac


def match_genes(
    genomad_genes: list[GenomadGene],
    bakta_features: dict[str, list[BaktaFeature]],
    contig_lengths: dict[str, int],
    min_overlap_frac: float = 0.5,
) -> list[dict]:
    """
    Attempt to match each geNomad gene to a Bakta CDS by coordinate overlap.

    For provirus genes, the geNomad coordinates are absolute positions within
    the parent contig. Bakta annotates the full genome, so contigs should match
    if the assembly is the same. We match on contig name substring.
    """
    results = []

    # Build a lookup: short contig key -> list of bakta features
    # Bakta contig names are like "contig_1", "contig_10", etc.
    # geNomad names may be:
    #   Bakta-style:  "21171_DNR38.contig_4.length_152347"  -> extract contig_4
    #   SPAdes-style: "NODE_10_length_67918_cov_257.546056" -> extract NODE number,
    #                 then map to contig_N by matching length in the GFF3

    import re

    def extract_contig_key(name: str, contig_map: dict = None) -> str:
        """
        Extract a Bakta contig key from various name formats.
        contig_map: optional dict mapping genomad contig strings -> bakta contig names
        """
        if contig_map and name in contig_map:
            return contig_map[name]
        # Standard contig_N format
        m = re.search(r"(contig_\d+)", name)
        if m:
            return m.group(1)
        return name

    # Build length -> bakta contig name mapping using exact ##sequence-region lengths
    length_to_bakta: dict[int, str] = {}
    for cname, clen in contig_lengths.items():
        if clen in length_to_bakta:
            print(f"  WARNING: duplicate contig length {clen} "
                  f"({length_to_bakta[clen]} and {cname}) - SPAdes auto-mapping may be ambiguous",
                  file=sys.stderr)
        length_to_bakta[clen] = cname

    # Auto-build contig map for SPAdes NODE names using embedded length
    def build_spades_map(genomad_genes_list, length_map) -> dict[str, str]:
        cmap = {}
        for gg in genomad_genes_list:
            name = gg.parent_contig
            if re.search(r"contig_\d+", name):
                continue  # already contig_N style
            # SPAdes: NODE_N_length_XXXX_cov_...
            m = re.search(r"NODE_\d+_length_(\d+)", name)
            if m:
                node_len = int(m.group(1))
                # Find Bakta contig with closest length
                # Exact match first
                if node_len in length_map:
                    cmap[name] = length_map[node_len]
                else:
                    # Closest length fallback (SPAdes length in name = full contig length,
                    # Bakta max coord may be slightly less due to annotation not reaching end)
                    closest = min(length_map.keys(), key=lambda x: abs(x - node_len))
                    if abs(closest - node_len) <= 50:  # within 50 bp tolerance
                        cmap[name] = length_map[closest]
                        print(f"  INFO: {name} (len {node_len}) -> {length_map[closest]} "
                              f"(len {closest}, diff {abs(closest-node_len)} bp)",
                              file=sys.stderr)
                    else:
                        print(f"  WARNING: could not auto-map {name} (len {node_len}) "
                              f"to any Bakta contig (closest: {closest})",
                              file=sys.stderr)
        return cmap

    auto_map = build_spades_map(genomad_genes, length_to_bakta)
    if auto_map:
        print(f"  SPAdes NODE -> Bakta contig auto-mappings:", file=sys.stderr)
        for k, v in auto_map.items():
            print(f"    {k} -> {v}", file=sys.stderr)

    bakta_by_contig: dict[str, list[BaktaFeature]] = {}
    for contig_name, feats in bakta_features.items():
        key = extract_contig_key(contig_name, auto_map)
        bakta_by_contig.setdefault(key, []).extend(feats)

    for gg in genomad_genes:
        contig_key = extract_contig_key(gg.parent_contig, auto_map)
        # Resolve through auto_map if SPAdes style
        if gg.parent_contig in auto_map:
            contig_key = auto_map[gg.parent_contig]

        # For provirus genes, coordinates are within the full contig
        # For standalone viral contigs, coords are within that short contig
        g_start = gg.start
        g_end = gg.end

        candidate_bakta = bakta_by_contig.get(contig_key, [])
        matched: list[BaktaFeature] = []

        for bf in candidate_bakta:
            if overlaps(g_start, g_end, bf.start, bf.end, min_overlap_frac):
                matched.append(bf)

        if matched:
            for bf in matched:
                results.append({
                    "genomad_gene_id": gg.gene_id,
                    "contig_key": contig_key,
                    "provirus_region": gg.provirus_region or "standalone",
                    "genomad_start": gg.start,
                    "genomad_end": gg.end,
                    "genomad_strand": gg.strand,
                    "genomad_marker": gg.marker,
                    "genomad_evalue": gg.evalue,
                    "genomad_bitscore": gg.bitscore,
                    "virus_hallmark": gg.virus_hallmark,
                    "genomad_description": gg.annotation_description,
                    "bakta_locus_tag": bf.locus_tag,
                    "bakta_product": bf.product,
                    "bakta_start": bf.start,
                    "bakta_end": bf.end,
                    "bakta_strand": bf.strand,
                    "match_status": "MATCHED",
                })
        else:
            results.append({
                "genomad_gene_id": gg.gene_id,
                "contig_key": contig_key,
                "provirus_region": gg.provirus_region or "standalone",
                "genomad_start": gg.start,
                "genomad_end": gg.end,
                "genomad_strand": gg.strand,
                "genomad_marker": gg.marker,
                "genomad_evalue": gg.evalue,
                "genomad_bitscore": gg.bitscore,
                "virus_hallmark": gg.virus_hallmark,
                "genomad_description": gg.annotation_description,
                "bakta_locus_tag": "NO_MATCH",
                "bakta_product": "NO_MATCH",
                "bakta_start": "NA",
                "bakta_end": "NA",
                "bakta_strand": "NA",
                "match_status": "UNMATCHED",
            })

    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gff3", required=True,
                        help="Bakta GFF3 annotation file")
    parser.add_argument("--genomad", required=True,
                        help="geNomad virus_genes.tsv file")
    parser.add_argument("--output", required=True,
                        help="Output TSV file")
    parser.add_argument("--min-overlap", type=float, default=0.5,
                        help="Minimum fractional overlap to call a match "
                             "(default: 0.5 = 50%% of shorter feature)")
    args = parser.parse_args()

    print(f"Parsing GFF3: {args.gff3}", file=sys.stderr)
    bakta, contig_lengths = parse_gff3(args.gff3)
    total_cds = sum(len(v) for v in bakta.values())
    print(f"  Loaded {total_cds} CDS features across {len(bakta)} contigs",
          file=sys.stderr)
    print(f"  Contigs in GFF3: {list(bakta.keys())}", file=sys.stderr)
    print(f"  Sequence-region lengths loaded for: {list(contig_lengths.keys())}",
          file=sys.stderr)

    print(f"Parsing geNomad: {args.genomad}", file=sys.stderr)
    genomad = parse_genomad(args.genomad)
    print(f"  Loaded {len(genomad)} viral genes", file=sys.stderr)

    print(f"Matching (min overlap: {args.min_overlap:.0%})...", file=sys.stderr)
    results = match_genes(genomad, bakta, contig_lengths, args.min_overlap)

    matched = sum(1 for r in results if r["match_status"] == "MATCHED")
    print(f"  {matched}/{len(results)} genes matched to Bakta CDS",
          file=sys.stderr)

    fieldnames = [
        "genomad_gene_id", "contig_key", "provirus_region",
        "genomad_start", "genomad_end", "genomad_strand",
        "genomad_marker", "genomad_evalue", "genomad_bitscore",
        "virus_hallmark", "genomad_description",
        "bakta_locus_tag", "bakta_product",
        "bakta_start", "bakta_end", "bakta_strand",
        "match_status",
    ]

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    print(f"Output written to: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
