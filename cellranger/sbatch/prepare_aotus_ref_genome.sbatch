#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=124GB
#SBATCH --job-name=prepare_aotus_ref_genome
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aet359@nyu.edu
#SBATCH --output=prepare_aotus_ref_genome_%j.out

# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

working_dir=$SLURM_SUBMIT_DIR

cd $working_dir

#STAR preparing genome using annotation file
cellranger mkref --genome=Aotus --fasta=aotus.fa --genes=aotus.filtered.gtf \
    --genome=Vivax --fasta=vivax.fa --genes=vivax.filtered.gtf

#--genome=will be the name of the output directory
#You can now specify this reference on the command line: cellranger --transcriptome=/scratch/aet359/scRNAseq/genomes/Aotus
cd ..

echo "Done, great job!"
