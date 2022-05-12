# DECIDEWP1_simulation

Repository for the simulation study side of WP1 on the DECIDE project. These scripts setup communities which are based on the butterfly and moth communities in the main project, to evaluate different adaptive sampling metrics. These include methods to approximate as well as possible those used in the main project. 

## Where's the code?

All code to run the simulations are found in `scripts/lotus/` used in the following order:

1. `slurm_simulate_species.R` 
2. `slurm_run_sim_sdm.R`
3. `slurm_adaptive_sample.R`
4. rerun `slurm_run_sim_sdm.R` with the AS methods added as parameters
5. `slurm_evaluate.R`

## Version explanation

- v1
- v2
- v3_asv1: community 1:50 run on 14/01/2022 with mechanism to sample across all species and with detectability set to 0.2. For adaptive sampling: `uptake = NULL`. Final filename `= asv1_v3combined_outputs_comm1_50_spp50.csv`
- v3_asv2: Same communities as v3, but adaptive sampling `uptake = 0.5`; filename `= asv2_v3combined_outputs_comm1_50_spp50.csv`. Started adaptive sampling 31/01/22, finished 1/2/22. 
- v3_asv3: Same communities as v3, but adaptive sampling `uptake = 1`. Started adaptive sampling 1/2/22, finished 2/2/2022. Filename `= asv3_v3combined_outputs_comm1_50_spp50.csv`.
- v3_asv4: Same communities as v3, but adaptive sampling `uptake = 0.1`. Started adaptive sampling running 2/2/22. Finished 3/2/2022. Filename `= asv4_v3combined_outputs_comm1_50_spp50.csv`.
- v3_asv5: Same communities as v3, but adaptive sampling `uptake = 0.01`. Started adaptive sampling running 4/2/22. Finished 9/2/22. Filename `= asv5_v3combined_outputs_comm1_50_spp50.csv`.
- v3_asv6: same communities as v3, models on communities 1:50. `uptake = 0.01`. Noticed mistake in adaptive sampling function that detectability is set to 0.2. This will effect results relative to 'none' method. Started AS (for all 50 communities) on 9/2/22, Finished 16.02.22
- v3_asv7: same communities as v3, models on communities 1:50. `uptake = 0.1`. Noticed mistake in adaptive sampling function that detectability is set to 0.2. This will effect results relative to 'none' method. Started AS (for all 50 communities) on 16.02.22
- v3_asv8: v3 communities, models on communities 1:50. `uptake = 0.1`. Detectability is set to 0.2. Removed `probability_weight_adj` parameter from the adaptive sampling methods with uptake as it was likely reinforcing biases.
- v4_asv1: v4 communities, models on communities 1:50 with fixed bugs to do with environmental variables being different between models. Detectability set to 0.2. Uptake set to `0.1`. Finished 12.5.22
- v4_asv2: uptake set to `0.01`. Started 12.5.22.
- v4_asv3: uptake set to `0`
- v4_asv4: uptake set to `0.5`



*To be continued...*
