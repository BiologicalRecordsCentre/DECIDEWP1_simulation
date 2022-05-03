
dirs <- config::get("LOTUSpaths")

# name of the community version we are running - so we're not overwriting things, keep same as for slurm_run_sim_sdm
community_version = 'v4'

# name of the adaptive sampling version we are looking to evaluate
AS_version = 'asv1'

# the name of the simulation run - same as slurm_simulate species
simulation_run_name = 'communities_1km'

n_communities = 1

n_species = 1:50

## index file
pars <- data.frame(community_folder = paste0(dirs$outpath, community_version, simulation_run_name, "/", community_version, sprintf("community_%i_%i_sim/", n_communities, max(n_species))),
                   model = "rf, gam, lr", 
                   method = "initial, none, uncertainty, prevalence, unc_plus_prev, unc_plus_recs, coverage",
                   community_version = community_version,
                   AS_version = AS_version)

community_folder = pars$community_folder[1]
model = pars$model[1]
method = pars$method[1]
community_version = pars$community_version[1]
AS_version = pars$AS_version[1]

source('scripts/slurm_extract_predictions_function.R')

## index file
ps <- extract_predictions(community_folder = paste0(dirs$outpath, community_version, simulation_run_name, "/", community_version, sprintf("community_%i_%i_sim/", n_communities, max(n_species))),
                          model = "rf, gam, lr", 
                          method = "initial, none, uncertainty, prevalence, unc_plus_prev, unc_plus_recs, coverage",
                          community_version = community_version,
                          AS_version = AS_version)
