#!/bin/bash
# ======================================================================
# Thesis Chapter 4, Table S4.1
# Maps Ptolemaea defence calls (PROKKA IDs) onto Bakta gene IDs.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/3_Transcriptomics/7_GenomeAnnotations/1_StrainGenomes/C67_bakta/get_ptolemaea_genes.sh
# ======================================================================

PROKKA_GFF="/mnt/scratch2/users/40204129/11_DefenseProfileofSsuis/6_FixingPtolemaea/Ptolemaea/output/01_prokka/C67_HBHC_prefixed/C67_HBHC_prefixed.gff"
BAKTA_GFF="/mnt/scratch2/users/40204129/3_Transcriptomics/7_GenomeAnnotations/1_StrainGenomes/C67_bakta/C67.gff3"
DEFENCE_CSV="/mnt/scratch2/users/40204129/11_DefenseProfileofSsuis/6_FixingPtolemaea/Ptolemaea/output/05_consensus/C67_HBHC_prefixed_defenceprofile.csv"

echo "=== Step 1: Extract defence gene IDs ==="
tail -n +2 "$DEFENCE_CSV" | cut -d',' -f1 | \
    sed 's/C67_HBHC_prefixed@//' > /tmp/defence_prokka_ids.txt
echo "Defence genes found: $(wc -l < /tmp/defence_prokka_ids.txt)"

echo ""
echo "=== Step 2: Extract PROKKA coordinates, normalise contig names ==="
# contig00001 -> contig_1, contig00012 -> contig_12 etc.
grep -F -f /tmp/defence_prokka_ids.txt "$PROKKA_GFF" | \
    grep -v "^#" | \
    awk '$3=="CDS"' | \
    awk '{
        # Normalise contig name: contig00001 -> contig_1
        contig = $1
        sub(/^contig0+/, "contig_", contig)
        match($9, /ID=([^;]+)/, id)
        print id[1] "\t" contig "\t" $4 "\t" $5 "\t" $7
    }' | sort -k2,2 -k3,3n > /tmp/defence_prokka_coords.tsv

echo "Coordinates extracted:"
cat /tmp/defence_prokka_coords.tsv

echo ""
echo "=== Step 3: Cross-reference with Bakta GFF ==="
echo -e "prokka_id\tcontig\tstart\tend\tstrand\tbakta_id\tbakta_product" > C67_defence_id_mapping.tsv

while IFS=$'\t' read -r prokka_id contig start end strand; do
    bakta_hit=$(awk -v s="$start" -v e="$end" -v c="$contig" '
        $3=="CDS" && $1==c &&
        int($4) >= int(s)-5 && int($4) <= int(s)+5 &&
        int($5) >= int(e)-5 && int($5) <= int(e)+5 {
            match($9, /ID=([^;]+)/, id)
            match($9, /Name=([^;]+)/, name)
            print id[1] "\t" name[1]
        }
    ' "$BAKTA_GFF")

    if [ -z "$bakta_hit" ]; then
        echo -e "$prokka_id\t$contig\t$start\t$end\t$strand\tNO_MATCH\tNO_MATCH" >> C67_defence_id_mapping.tsv
    else
        echo -e "$prokka_id\t$contig\t$start\t$end\t$strand\t$bakta_hit" >> C67_defence_id_mapping.tsv
    fi
done < /tmp/defence_prokka_coords.tsv

echo ""
echo "=== Final mapping ==="
column -t C67_defence_id_mapping.tsv

# Quick sanity check
matched=$(grep -v "NO_MATCH" C67_defence_id_mapping.tsv | tail -n +2 | wc -l)
total=$(tail -n +2 C67_defence_id_mapping.tsv | wc -l)
echo ""
echo "Matched: $matched / $total"
