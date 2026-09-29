#!/usr/bin/env python3
# ======================================================================
# Thesis Chapter 4, Methods: genome annotation
# Merges Bakta, eggNOG, geNomad and Ptolemaea calls into one GFF3 per strain.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/3_Transcriptomics/7_GenomeAnnotations/3_ConsolidatingAnnotations/merge_annotations.py
# ======================================================================
"""
Merge multiple annotation sources into a single GFF3 file
"""

import sys
import csv
from collections import defaultdict

def load_defensefinder(tsv_file):
    """Load DefenseFinder annotations"""
    annotations = {}
    with open(tsv_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            locus = row['hit_id']
            # Format: DefenseFinder:type|subtype|activity
            annot = f"DefenseFinder:{row['type']}|{row['subtype']}|{row['activity']}"
            annotations[locus] = annot
    return annotations

def load_padloc(csv_file):
    """Load PADLOC annotations"""
    annotations = {}
    with open(csv_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            locus = row['target.name']
            # Format: PADLOC:system
            annot = f"PADLOC:{row['system']}"
            annotations[locus] = annot
    return annotations

def load_eggnog(emapper_file):
    """Load eggNOG functional annotations"""
    annotations = {}
    skipped = 0
    total_lines = 0
    
    with open(emapper_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            total_lines += 1
            parts = line.strip().split('\t')
            
            if len(parts) < 20:
                skipped += 1
                print(f"  Skipped line {total_lines}: only {len(parts)} columns")
                continue
            
            locus = parts[0]
            description = parts[7] if parts[7] != '-' else ''
            preferred_name = parts[8] if parts[8] != '-' else ''
            cog_category = parts[6] if parts[6] != '-' else ''
            kegg_ko = parts[11] if parts[11] != '-' else ''
            kegg_pathway = parts[12] if parts[12] != '-' else ''
            ec = parts[10] if parts[10] != '-' else ''
            
            # Build eggNOG annotation string
            annot_parts = []
            if preferred_name:
                annot_parts.append(f"gene:{preferred_name}")
            if cog_category:
                annot_parts.append(f"COG:{cog_category}")
            if ec:
                annot_parts.append(f"EC:{ec}")
            if kegg_ko:
                annot_parts.append(f"KO:{kegg_ko}")
            if kegg_pathway:
                pathways = kegg_pathway.split(',')[:3]
                annot_parts.append(f"pathway:{','.join(pathways)}")
            if description:
                desc = description[:100] + '...' if len(description) > 100 else description
                annot_parts.append(f"function:{desc}")
            
            if annot_parts:
                annotations[locus] = "eggNOG:" + "|".join(annot_parts)
    
    print(f"  Total eggNOG lines processed: {total_lines}")
    print(f"  Lines skipped (< 20 columns): {skipped}")
    
    return annotations

def parse_gff_attributes(attr_string):
    """Parse GFF3 attributes into dictionary"""
    attrs = {}
    for item in attr_string.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            attrs[key] = value
    return attrs

def format_gff_attributes(attrs):
    """Format attributes dictionary back to GFF3 string"""
    # Maintain order: ID, Name, locus_tag, product, Dbxref, Note
    ordered_keys = ['ID', 'Name', 'locus_tag', 'product', 'transl_table', 'Dbxref', 'Note']
    result = []
    
    # Add ordered keys first
    for key in ordered_keys:
        if key in attrs:
            result.append(f"{key}={attrs[key]}")
    
    # Add remaining keys
    for key, value in attrs.items():
        if key not in ordered_keys:
            result.append(f"{key}={value}")
    
    return ';'.join(result)

def merge_gff(input_gff, output_gff, *annotation_dicts):
    """Merge all annotations into GFF3"""
    
    stats = {'total': 0, 'annotated': 0}
    
    with open(input_gff) as inf, open(output_gff, 'w') as outf:
        for line in inf:
            # Pass through header lines
            if line.startswith('#'):
                outf.write(line)
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                outf.write(line)
                continue
            
            # Only process CDS features
            if fields[2] != 'CDS':
                outf.write(line)
                continue
            
            stats['total'] += 1
            
            # Parse attributes
            attrs = parse_gff_attributes(fields[8])
            locus_tag = attrs.get('locus_tag', '')
            
            if not locus_tag:
                outf.write(line)
                continue
            
            # Collect annotations from all sources
            notes = []
            for annot_dict in annotation_dicts:
                if locus_tag in annot_dict:
                    notes.append(annot_dict[locus_tag])
            
            # Add annotations to Note field
            if notes:
                stats['annotated'] += 1
                if 'Note' in attrs:
                    # Bakta might already have Note field
                    attrs['Note'] = attrs['Note'] + '|' + '|'.join(notes)
                else:
                    attrs['Note'] = '|'.join(notes)
            
            # Write updated line
            fields[8] = format_gff_attributes(attrs)
            outf.write('\t'.join(fields) + '\n')
    
    return stats

if __name__ == '__main__':
    if len(sys.argv) != 6:
        print("Usage: merge_annotations.py input.gff3 defensefinder.tsv padloc.csv eggnog.emapper output.gff3")
        sys.exit(1)
    
    input_gff = sys.argv[1]
    df_file = sys.argv[2]
    padloc_file = sys.argv[3]
    eggnog_file = sys.argv[4]
    output_gff = sys.argv[5]
    
    print("Loading annotations...")
    df_annots = load_defensefinder(df_file)
    print(f"  DefenseFinder: {len(df_annots)} annotations")
    
    padloc_annots = load_padloc(padloc_file)
    print(f"  PADLOC: {len(padloc_annots)} annotations")
    
    eggnog_annots = load_eggnog(eggnog_file)
    print(f"  eggNOG: {len(eggnog_annots)} annotations")
    
    print("\nMerging annotations into GFF3...")
    stats = merge_gff(input_gff, output_gff, df_annots, padloc_annots, eggnog_annots)
    
    print(f"\nDone!")
    print(f"  Total CDS features: {stats['total']}")
    print(f"  Features with new annotations: {stats['annotated']}")
    print(f"  Output written to {output_gff}")

