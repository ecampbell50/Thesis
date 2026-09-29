#!/bin/bash
# ======================================================================
# Thesis Chapter 5, Methods: prophage annotation
# geNomad end-to-end (--conservative --sensitivity 5) on all 2,119 genomes.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/Chapter3_Coevolution/RunGeNomad_24Mar26.sub
# ======================================================================

#SBATCH --job-name=geNomad_24May26
#SBATCH --error=geNomad_24May26.err
#SBATCH --output=geNomad_24May26.out
#SBATCH --partition=k2-bioinf,k2-lowpri
#SBATCH --time=2-23:59:59
#SBATCH --nodes=1
#SBATCH --mem=128G
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
