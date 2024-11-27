#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G
#SBATCH --time=22:00:00
#SBATCH --mail-user=cdavilaojed18@students.kgi.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="md_gpu_final"

# Print current date
date

# Load gromacs
module load gromacs/2024.3-gpu

#run the application
gmx mdrun -deffnm md_0_1 -nstlist 40 -nb gpu
