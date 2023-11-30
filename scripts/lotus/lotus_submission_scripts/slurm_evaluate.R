
library(rslurm)

source("scripts/slurm_evaluate_function.R")

dirs <- config::get("LOTUSpaths")

# name of the community version we are running - so we're not overwriting things, keep same as for slurm_run_sim_sdm
community_version = 'v4'

as_versions = c('asv1', 'asv2', 'asv3', 'asv4')

for(i in 1:4) {
  
  # name of the adaptive sampling version we are looking to evaluate
  AS_version = as_versions[i]
  
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
  
  #### slurm apply call
  sdm_slurm <- slurm_apply(slurm_evaluate,
                           params = pars,
                           jobname = paste('evaluate_community', community_version, AS_version, collapse = '_'),
                           nodes = length(pars$community_folder),
                           cpus_per_node = 1,
                           slurm_options = list(partition = 'short-serial-4hr',
                                                time = '3:59:59',
                                                mem = 3000,
                                                output = "sim_eval_%a.out",
                                                error = "sim_eval_%a.err",
                                                account = "short4hr"),
                           sh_template = "jasmin_submit_sh.txt",
                           submit = T)
  
}