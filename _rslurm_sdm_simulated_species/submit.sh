#!/bin/bash
#
#SBATCH --array=0-9
#SBATCH --cpus-per-task=1
#SBATCH --job-name=sdm_simulated_species
#SBATCH --output=slurm_%a.out
#SBATCH --partition=short-serial
#SBATCH --time=23:59:59
#SBATCH --mem=20000
module add jaspy
 Rscript slurm_run.R
