source("slurm_run_sim_sdm.R")

index <- as.numeric(commandArgs(trailingOnly = TRUE))[1]

slurm_run_sim_sdm(index)
