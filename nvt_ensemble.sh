#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --time=01:00:00
#SBATCH --mail-user=cdavilaojed18@students.kgi.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Ensemble 1 NVT"

# Print current date
date

# Load gromacs
module load gromacs

#run the application
gmx mdrun -deffnm nvt
