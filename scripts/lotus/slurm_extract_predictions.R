
dirs <- config::get("LOTUSpaths")
library(rslurm)

# name of the community version we are running - so we're not overwriting things, keep same as for slurm_run_sim_sdm
community_version = 'v4'

# name of the adaptive sampling version we are looking to evaluate
AS_version = 'asv2'

# the name of the simulation run - same as slurm_simulate species
simulation_run_name = 'communities_1km'

n_communities = 1:50

n_species = 1:50

## index file
pars <- data.frame(community_folder = paste0(dirs$outpath, community_version, simulation_run_name, "/", community_version, sprintf("community_%i_%i_sim/", n_communities, max(n_species))),
                   model = "rf, gam, lr", 
                   method = "initial, none, uncertainty, prevalence, unc_plus_prev, unc_plus_recs, coverage",
                   community_version = community_version,
                   AS_version = AS_version)

source('scripts/slurm_extract_predictions_function.R')

# ## index file
# ps <- extract_predictions(community_folder = paste0(dirs$outpath, community_version, simulation_run_name, "/", community_version, sprintf("community_%i_%i_sim/", n_communities, max(n_species)))[1],
#                           model = "rf, gam, lr", 
#                           method = "initial, none, uncertainty, prevalence, unc_plus_prev, unc_plus_recs, coverage",
#                           community_version = community_version,
#                           AS_version = AS_version)

#### slurm apply call
sdm_slurm <- slurm_apply(extract_predictions,
                         params = pars,
                         jobname = paste('predictions_community', community_version, AS_version, collapse = '_'),
                         nodes = length(pars$community_folder),
                         cpus_per_node = 1,
                         slurm_options = list(partition = 'short-serial',
                                              time = '23:59:59',
                                              mem = 3000,
                                              output = "sim_eval_%a.out",
                                              error = "sim_eval_%a.err"),
                         sh_template = "jasmin_submit_sh.txt",
                         submit = T)



for(i in n_communities) {
  
  print(i)
  ps <- extract_predictions(community_folder = pars$community_folder[i],
                            model = pars$model[i], 
                            method = pars$method[i],
                            community_version = pars$community_version[i],
                            AS_version = pars$AS_version[i])
  
}
