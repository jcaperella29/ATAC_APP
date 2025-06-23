#!/bin/bash
#SBATCH --job-name=ATAC_Shiny
#SBATCH --output=atac_shiny.log
#SBATCH --error=atac_shiny.err
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=interactive

# Load Singularity module (may vary by cluster)
module load singularity

# Set working directory
cd $SLURM_SUBMIT_DIR

# Define the Singularity image and bind paths
SIF="atac-shiny.sif"
APP_DIR=$(pwd)

# Run the container
singularity run --bind $APP_DIR:/mnt $SIF
