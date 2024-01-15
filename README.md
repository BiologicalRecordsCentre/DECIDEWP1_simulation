# Simulating adaptive sampling by citizen scientists  

This repository contains the R code needed to recreate the analyses in the paper "Adaptive sampling by citizen scientists improves species distribution model performance: a simulation study", part of the NERC-funded DECIDE project. In these simulations, we investigated the effect of different methods of adaptive sampling on the perfomance of species distribution models. We also explored the effect of simulated recorder participation (what we call uptake) on the effect of adaptive sampling on model performance.  

This code is written to be run on a SLURM cluster to enable easy modelling of multiple species of interest and relies heavily on the rslurm package to submit jobs (Marchaud _et al._ 2023 https://cran.r-project.org/web/packages/rslurm/index.html). The necessary scripts need to be transfered to the chosen cluster and the locations of these need to be specified using config files (Allaire _et al._ 2023 https://cran.r-project.org/web/packages/config/index.html). 

## Where's the code?

All code can be found in the `scripts/` folder.

```bash

|- scripts
|  ├─ config_files
|  ├─ evaluation
|  ├─ figures
│  ├─ lotus
|    ├ lotus_functions
|    ├ lotus_other
|    ├ lotus_submission_scripts

```

The scripts to run the main analyses on a SLURM cluster can be found in the `lotus/` folder. The functions used to carry out the various parts of the analysis are in the `functions/` folder. Code to plot the figures in the main manuscript and supplementary materials can be found in `figures/` folder. The `evaluation` folder contains scripts to perform various things associated with checking the main results and several aspects of the simulation workflow. These are not used to perform any of the analysis or plotting in the paper. The `scripts/processing_environmental_data.R` script was used to process all of the environmental layers into a format useable in the analyses. The processed environmental data used to run all the analyses are provided. 

### The main analysis

The code to run the main analyses on the SLURM cluster are split into three folders. The scripts to submit jobs to the SLURM cluster can all be found in the folder `lotus/lotus_submission_scripts/`. Each of the scripts in here corresponds to one of the scripts in the folder `lotus/lotus_functions/`, which contains the functions used by `rslurm::slurm_apply()` to run the various parts of the analyses. Once all scripts are in their correct locations, scripts in `lotus/lotus_submission_scripts/` need to be run in the following order to carry out the adaptive sampling workflow.

1. `slurm_simulate_species.R`

Simulate `n` species in each of the `n_communities`. Outputs communities of `n` species as Rdata files.


2. `slurm_run_sim_sdm.R`

Run the species distribution models. Set the `data_type` parameter to `initial_AS_coverage` for the first model run to get the baseline models.


3. `slurm_adaptive_sample.R`

Carry out the adaptive sampling using the methods described in the paper.


4. rerun `slurm_run_sim_sdm.R`

Run the species distribution models for a second time. This time setting the  `data_type` parameter to a vector containing the names of the various adaptive sampling methods.


5. `slurm_evaluate.R`

Carry out the evaluation comparing the baseline and adaptively sampled species distribution models to the true distributions of each species.

This completes the main analyses which are run on the SLURM cluster. The folder `lotus_other` contains scripts to perform various other aspects useful for the paper. Most importantly, the script `combine_files_on_lotus.R` combines all of the evaluation outputs into a single file for subsequent analysis.
