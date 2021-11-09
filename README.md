# DECIDEWP1_simulation

Repository for the simulation study side of WP1 on the DECIDE project. These scripts setup communities which are based on the butterfly and moth communities in the main project, to evaluate different adaptive sampling metrics. These include methods to approximate as well as possible those used in the main project. 

## Where's the code?

All code to run the simulations are found in `scripts/lotus/` used in the following order:

1. `slurm_simulate_species.R` 
2. `slurm_run_sim_sdm.R`
3. `slurm_adaptive_sample.R`
4. rerun `slurm_run_sim_sdm.R` with the AS methods added as parameters
5. `slurm_evaluate.R`

*To be continued...*
