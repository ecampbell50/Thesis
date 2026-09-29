#!/bin/bash
# ======================================================================
# Thesis Chapter 3, Methods: Serotype (Table S4)
# Runs SsuisSerotyping_pipeline on the pseudo-FASTQs (BBMap conversion step was run by hand, not scripted).
# Original location (Kelvin2 HPC): /mnt/autofs/mcclayrds-projects/r3298gfs/2_After29Nov23/code/1_sbatches/AllGenomesSerotyping.sub
# ======================================================================

#SBATCH --job-name=AllGenomesSerotyping
#SBATCH --error=AllGenomesSerotyping.err
#SBATCH --output=AllGenomesSerotyping.log
#SBATCH --partition=k2-lowpri,bio-compute
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

perl Ssuis_serotypingPipeline.pl --fastq_directory /users/40204129/sharedscratch/13_Serotyping/10_AllGenomesRun1/fastqs --scoreName AllGenomes --ends pe
