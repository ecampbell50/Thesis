#!/bin/bash
# ======================================================================
# Thesis Chapter 5, Methods: prophage annotation
# Writes one FASTA per provirus (input for sourmash and the spacer BLAST).
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/Chapter3_Coevolution/04_genomad/ExtractProviruses.sh
# ======================================================================

out_dir="proviruses"
mkdir -p "$out_dir"

echo "Extracting proviruses..."

for genomad_dir in *_genomad; do
    # Extract genome ID directly from directory name (most reliable)
    genome_id="${genomad_dir%_genomad}"
    
    virus_file="${genomad_dir}/${genome_id}_summary/${genome_id}_virus.fna"
    
    if [[ ! -f "$virus_file" ]]; then
        continue
    fi
    
    echo "Processing $genome_id..."
    
    # Use awk with RS=">" to split on FASTA headers
    awk -v gid="$genome_id" -v outdir="$out_dir" '
    BEGIN { RS=">" }
    NF && $1 ~ /\|provirus_/ {
        header = $1
        
        # Remove leading > if present
        gsub(/^>/, "", header)
        
        # Parse: CONTIG_ID|provirus_START_END
        split(header, parts, "|")
        contig = parts[1]
        provirus_info = parts[2]
        
        # Replace underscores with dashes in coordinates
        gsub(/_/, "-", provirus_info)
        
        # Build filename
        filename = outdir "/" gid "_" contig "_" provirus_info ".fna"
        
        # Write header and sequence
        print ">" header > filename
        for (i=2; i<=NF; i++) {
            print $i >> filename
        }
    }
    ' "$virus_file"
done

# Fix sequence headers to include genome ID
for file in proviruses/*.fna; do
    filename=$(basename "$file")
    genome_id=$(echo "$filename" | cut -d_ -f1)
    contig=$(echo "$filename" | cut -d_ -f2)
    
    sed -i "s/^>.*|provirus_/>$genome_id\_$contig|provirus_/" "$file"
    sed -i "s/|/_/g" "$file"
done


echo ""
echo "✅ Done!"
echo "Sample files created:"
ls -1 "$out_dir" | head -5
echo "... ($(ls -1 "$out_dir" | wc -l) total files)"
