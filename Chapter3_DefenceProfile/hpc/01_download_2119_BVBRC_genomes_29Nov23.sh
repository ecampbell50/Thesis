#!/bin/bash
# ======================================================================
# Thesis Chapter 3, Methods: Genome collection
# Downloads the 2,119 S. suis assemblies from the BV-BRC FTP (29 Nov 2023).
# Original location (Kelvin2 HPC): /mnt/autofs/mcclayrds-projects/r3298gfs/2_After29Nov23/code/1_sbatches/wgetGenomes_29Nov23.sub
# ======================================================================

#SBATCH --job-name=wgetGenomes_29Nov23
#SBATCH --error=wgetGenomes_29Nov23.error
#SBATCH --output=wgetGenomes_29Nov23.log
#SBATCH --partition=bio-compute,k2-lowpri
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

echo "Today is Wednesday 29 Nov 2023..."
echo "Starting download of genomes..."

input_file="/users/40204129/sharedscratch/1_PangenomeAnalysis/2_After29Nov23/data/Ssuis_GenomeIDs_29Nov23.csv"
output_folder="/users/40204129/sharedscratch/1_PangenomeAnalysis/2_After29Nov23/data/1_GenomesFromBVBRC"

# Loop through each genome ID in the input file
while IFS=, read -r genome_id; do
    # Remove double quotes from the genome ID
    genome_id=${genome_id//\"/}

    # Download the genome from the FTP site
    wget -qN "ftp://ftp.bvbrc.org/genomes/$genome_id/$genome_id.fna" -P "$output_folder"

    # Check if the genome was downloaded
    if [ $? -eq 0 ]; then
        echo "$genome_id successfully downloaded..."
    else
        echo "Error: $genome_id not downloaded..."
    fi
done < "$input_file"
