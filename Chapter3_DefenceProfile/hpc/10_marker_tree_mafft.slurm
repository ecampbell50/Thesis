#!/bin/bash
# ======================================================================
# Thesis Chapter 3, Figure 3.7 (marker-gene species tree)
# MAFFT alignment of each COG.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/11_DefenseProfileofSsuis/1_SsuisGenomes/ssuis_defence_project/GeneReconcilliation_of_tmn/COG_to_Seq_BLASTS/MAFFT_alignment.sub
# ======================================================================
#SBATCH --job-name=MAFFT_COGs
#SBATCH --error=MAFFT_COGs_%A_%a.err
#SBATCH --output=MAFFT_COGs_%A_%a.out
#SBATCH --partition=k2-sandbox
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --ntasks=2
#SBATCH --array=1-40
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Create output directory
mkdir -p alignments

# Get array of input files
FILES=(COG*_seqs.fas)

# Get the file for this array task (index starts at 0, array starts at 1)
INPUT_FILE="${FILES[$SLURM_ARRAY_TASK_ID-1]}"

# Extract COG name
COG_NAME=$(basename "$INPUT_FILE" _seqs.fas)

echo "Processing ${COG_NAME}..."
echo "Input file: ${INPUT_FILE}"
echo "Started at: $(date)"

# Run MAFFT with high accuracy for phylogenetics
mafft --maxiterate 1000 --localpair --thread 2 "$INPUT_FILE" > "alignments/${COG_NAME}_aligned.faa"

echo "Finished at: $(date)"
echo "Output: alignments/${COG_NAME}_aligned.faa"

