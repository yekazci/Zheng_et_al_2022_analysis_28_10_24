#!/bin/bash

# Submit this script inside parent folder of python_scripts.

# SLURM directives:

#SBATCH -J scvi-integrate
#SBATCH --export=ALL
#SBATCH --cpus-per-task=64
#SBATCH --mem=512G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=yusufenes.kazci@mdc-berlin.de
#SBATCH --output=slurm_output_%j.txt
#SBATCH --error=slurm_error_%j.txt
#SBATCH --chdir=./                # Set the working directory to the current, which is default.


srun python python_scripts/scvi_batch_integration.py

echo "completed."