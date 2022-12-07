#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=124GB
#SBATCH --job-name=cellranger_count_74
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aet359@nyu.edu
#SBATCH --output=cellranger_count_74_%j.out

# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

working_dir=$SLURM_SUBMIT_DIR

cd $working_dir

# Making a cell ranger count for each individual

module load cellranger/7.0.0
module load bcl2fastq/2.20.0.422

cellranger count --id=Aotus_cell_ranger_74  \
    --sample=SRR11008274 \
    --fastqs=data/aotus/SRR11008274_S3_R1/ \
    --transcriptome=/scratch/aet359/scRNAseq/genomes/Aotus_and_Vivax \
    --include-introns=true \
    --chemistry=SC5P-R2
    
#--id=made up name for the experiment



