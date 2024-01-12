# Simulating adaptive sampling by citizen scientists  

This repository contains the R code needed to recreate the analyses in the paper "Adaptive sampling by citizen scientists improves species distribution model performance: a simulation study", part of the NERC-funded DECIDE project. In these simulations, we investigated the effect of different methods of adaptive sampling on the perfomance of species distribution models. We also explored the effect of simulated recorder participation (what we call uptake) on the effect of adaptive sampling on model performance.  

This code is written to be run on a SLURM cluster to enable easy modelling of multiple species of interest and relies heavily on the rslurm package to submit jobs (Marchaud _et al._ 2023 https://cran.r-project.org/web/packages/rslurm/index.html). The necessary scripts need to be transfered to the chosen cluster and the locations of these need to be specified using config files (Allaire _et al._ 2023 https://cran.r-project.org/web/packages/config/index.html). 

## Where's the code?

Here is a brief overview of the main scripts and folder structure. 

All code to run the simulations are found in `scripts/lotus/` used in the following order:

1. `slurm_simulate_species.R` 
2. `slurm_run_sim_sdm.R`
3. `slurm_adaptive_sample.R`
4. rerun `slurm_run_sim_sdm.R` with the AS methods added as parameters
5. `slurm_evaluate.R`

The code to produce all figures can be found in `scripts/paper_scripts/figures/`.




*To be continued...*
