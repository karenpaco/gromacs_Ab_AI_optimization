#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=01:00:00
#SBATCH --mail-user=cdavilaojed18@students.kgi.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Energy Minimization"

# Print current date
date

# Load gromacs
module load gromacs

#run the application
gmx mdrun -v -deffnm em
