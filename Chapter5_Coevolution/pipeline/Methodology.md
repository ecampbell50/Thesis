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

***Steps 5-6***
Run geNomad and extract the proviruses: `hpc/01`-`03`
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

*Section 2 - CRISPR spacers*
MinCED and the spacer BLAST: `hpc/05`-`07`
