#!/bin/bash
#SBATCH --job-name=atac_app
#SBATCH --output=atac_app_%j.log
#SBATCH --error=atac_app_%j.err
#SBATCH --partition=debug
#SBATCH --mem=14000
#SBATCH --time=02:00:00

# Set working directory
cd $SLURM_SUBMIT_DIR

# Define the Singularity image and bind paths
SIF="atac-shiny.sif"
APP_DIR=$(pwd)

# Define the Shiny port to expose (match this with browser access)
export PORT=8787

# Debug print to confirm launch
echo "📣 Launching Shiny app on port $PORT..."
echo "📁 App dir: $APP_DIR"
echo "📦 Container: $SIF"

# Run the container with port env set
singularity run --env PORT=8787 --bind $APP_DIR:/app $SIF

# If you want to leave a message when it exits
echo "✅ ATAC Shiny container exited at $(date)"
