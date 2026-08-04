*Section 1 - Getting Bipartite Network of Bacteria and Prophages*

***Step 1:***
Retrieve 2,119 S. suis genomes from Kelvin storage
```bash
cp /mnt/autofs/mcclayrds-projects/r3298gfs/Ptolemaea/data/raw/001_GenomesfromBV-BRC/*fna /mnt/scratch2/users/40204129/Chapter3_Coevolution/genomes
```

***Step 2***
Run sourmash on genomes
sourmash (doi: 10.21105/joss.06830)

```bash
# Create new sourmash env cos pip ins't working anymore...
mamba create -n Sourmash_23May26 -c bioconda -c conda-forge sourmash
mamba activate Sourmash_23May26
sourmash --version

sourmash 4.9.4
```
***Step 3***
Run sourmash on genomes

```bash
srun --pty --partition=k2-bioinf --time=06:00:00 --cpus-per-task=4 --ntasks=1 --mem-per-cpu=16G bash

# Create sourmash sketches
# Use kmer size 21 to get some sensitivity
sourmash sketch dna -p scaled=1000,k=21 *fna
sourmash compare --ksize 21 --csv Ssuis_2119_k21s1000.csv *sig
```

***Step 4***
Run netpass yay

```bash
git clone https://github.com/ecampbell50/Netpass.git

mamba env create --name Netpass_23May26 --file Netpass/Netpass_env.yaml
2
conda activate Netpass_23May26

python3 Netpass/UpdatedNetpass.py network -i Ssuis_2119_k21s1000.csv -o results/ -t BAC
```

***Step 5***
Run geNomad

```bash
mamba create -n genomad_23may26 -c conda-forge -c bioconda genomad

mamba activate genomad_23May26

genomad download-database .
[22:20:04] geNomad database (v1.9) is ready to be use

#!/bin/bash

#SBATCH --job-name=geNomad_23May26
#SBATCH --error=geNomad_23May26.err
#SBATCH --output=geNomad_23May26.out
#SBATCH --partition=k2-bioinf,k2-lowpri
#SBATCH --time=2-23:59:59
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --ntasks=24
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

mkdir -p 04_genomad

for i in 01_genomes/*fna; do
        base=$(basename ${i%.fna})
        genomad end-to-end --conservative --restart --threads 24 --sensitivity 5 \
        $i 04_genomad/${base}_genomad genomad_db

        echo "Finished $i"
done
```
***Step 6***
Extract provirus sequences

```bash
#!/bin/bash

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
```
3,216 proviruses extracted

***Step 7***
Move proviruses to folder for sourmash

Mkdir 05_proviruses
mv 04_genomad/proviruses/* 05_proviruses

mkdir 06_provirus_sketches

```bash
srun --pty --partition=k2-bioinf,k2-hipri --time=03:00:00 --cpus-per-task=4 --ntasks=1 --mem-per-cpu=16G bash

# Kmer size 21 good for sensitivity/specificity of bacterial genomes
# Kmer size 15 for phages = more sensitive for divergent sequences
# Scaled = 500 = more sensitive/ larger sketches
sourmash sketch dna -p scaled=500,k=15 *fna
sourmash compare --ksize 15 --csv Ssuis_3216proviruses_k15s500.csv *sig
```

***Step 8***
Run netpass on proviruses

mkdir 07_provirus_netpass

```bash
python3 Netpass/UpdatedNetpass.py network -i Ssuis_3216proviruses_k15s500.csv -o results/ -t PRO
```

***Step 9***
Create genome-level bipartite network

```bash
# Link genome ID to provirus ID
# This outputs Source,Target,Value format to append the three edgetables together
for i in $(ls 05_proviruses/); do PROVIRUS=$i; GENOME=$(echo ${i/_*.fna/}); echo "${GENOME}.fna,$PROVIRUS,1"; done > 08_bipartite_network/Genome_to_Prophage_link.csv

# Create a genome-level bipartite network of bacteria and their prophages
        # Removes self-loops
        cat 03_netpass/results/Edgetable_Variations/Ssuis_2119_k21s1000_49.csv | awk -F',' '$1 != $2' >        08_bipartite_network/BAC_genomes_49_edgetable.csv
        # Remove self-loops and the header (will only need one header)
        cat 07_provirus_netpass/results/Edgetable_Variations/Ssuis_3216proviruses_k15s500_44.csv | awk -F',' '$1 != $2 && $1 != "Source"' > 08_bipartite_network/PRO_genomes_44_edgetable.csv

cat BAC_genomes_49_edgetable.csv > Bipartite_BAC_PRO_genome_edgetable.csv
cat PRO_genomes_44_edgetable.csv >> Bipartite_BAC_PRO_genome_edgetable.csv
cat Genome_to_Prophage_link.csv >> Bipartite_BAC_PRO_genome_edgetable.csv

# Get node table for cytoscape
for i in $(ls 01_genomes/); do echo "$i,BAC"; done > 08_bipartite_network/BAC_nodetable.csv
for i in $(ls 05_proviruses/); do echo "$i,PRO"; done > 08_bipartite_network/PRO_nodetable.csv
```
**MILESTONE 1: Bipartite Network Figure**

*Section 2 - Getting CRISPR spacers of all genomes to check past immunity*

Using minced v0.4.2 which is in the prokka env in igfs-anaconda

```bash
#SBATCH --job-name=crispr_spacers
#SBATCH --partition=k2-bioinf,k2-medpri
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=16G
#SBATCH --output=minced_17jun26.out
#SBATCH --error=minced_17jun26.err
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=END,FAIL

cd /mnt/scratch2/users/40204129/Chapter3_Coevolution/09_minced
: > all_spacers.fasta
for f in ../01_genomes/*.fna; do
  gid=$(basename "$f" .fna)
  minced -gffFull -spacers "$f" >/dev/null 2>&1      # spacers -> ./<gid>_spacers.fa ; text table dumped
  sp="${gid}_spacers.fa"
  [[ -s "$sp" ]] && sed "s/^>/>${gid}|/" "$sp" >> all_spacers.fasta
done
echo "spacers: $(grep -c '^>' all_spacers.fasta) from $(grep '^>' all_spacers.fasta | cut -d'|' -f1 | sort -u | wc -l) genomes"
```

"all_spacers" didn't work, and they don't have genome ID

```bash
cd /mnt/scratch2/users/40204129/Chapter3_Coevolution/09_minced
: > all_spacers.fasta
for sp in arrays/*_spacers.fa; do
  [[ -s "$sp" ]] || continue                 # skip empties (genomes with no arrays)
  gid=$(basename "$sp" _spacers.fa)          # e.g. 1005042.3
  sed "s/^>/>${gid}|/" "$sp" >> all_spacers.fasta
done
# sanity check:
echo "spacers: $(grep -c '^>' all_spacers.fasta) from $(grep '^>' all_spacers.fasta | cut -d'|' -f1 | sort -u | wc -l) genomes"
head -4 all_spacers.fasta
```

spacer prophage BLAST

```bash
#!/bin/bash
#SBATCH --job-name=spacer_blast
#SBATCH --partition=k2-bioinf,k2-lowpri
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=16G
#SBATCH --output=spacer_blast.out
#SBATCH --error=spacer_blast.err
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=END,FAIL

cd /mnt/scratch2/users/40204129/Chapter3_Coevolution/09_minced
BASE=/mnt/scratch2/users/40204129/Chapter3_Coevolution
# provirus db, with each provirus's FILENAME as its blast id (so hits map to prophage_id)
: > all_proviruses.fna
for f in "$BASE"/05_proviruses/*.fna; do
  id=$(basename "$f")
  awk -v id="$id" '/^>/{print ">"id; next}{print}' "$f" >> all_proviruses.fna
done
makeblastdb -in all_proviruses.fna -dbtype nucl -out provirus_db >/dev/null

# match spacers -> proviruses, then keep only protospacer-grade hits
blastn -task blastn-short -query all_spacers.fasta -db provirus_db -dust no -evalue 1 -num_threads 8 \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend sstart send evalue bitscore' \
  > spacer_vs_provirus.raw.tsv
awk -F'\t' '($3>=95)&&($5<=1)&&($4/$7>=0.90)' spacer_vs_provirus.raw.tsv > spacer_vs_provirus.tsv
echo "protospacer-grade hits: $(wc -l < spacer_vs_provirus.tsv)"
```

