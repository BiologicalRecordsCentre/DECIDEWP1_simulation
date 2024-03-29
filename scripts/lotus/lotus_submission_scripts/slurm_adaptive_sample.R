## script to submit adaptive sampling

# source function to run
source("scripts/slurm_adaptive_sample_function.R")

library(rslurm)

dirs <- config::get("LOTUSpaths_AS")

# name of the versions we are running - so we're not overwriting things
# Three different versions, one for community-level which includes the community folders and species models folders
community_version = 'v4'

## run multiple adaptive sampling versions at once

asv_version = data.frame(as_version =  c("asv1", "asv2", "asv3", "asv4"),
                         uptake_value = c(0.1, 0.01, 0, 0.5))

for(asv in 1:nrow(asv_version)){
  
  # and an adaptive sampling version, which is if we want to run the adaptive sampling 
  # process more than once - these outputs are stored in the same place as the old outputs
  # must always be prefixed by asv
  AS_version = asv_version$as_version[asv]
  
  # the name of the simulation run - same as slurm_simulate species
  simulation_run_name = 'communities_1km'
  
  # number of communities - a vector!
  n_communities = 1:50
  
  # number of species in each community - used only in the parameter file to allow runs with different numbers of species
  n_species = 1:50
  
  # the adaptive sampling methods to use 
  method = "coverage" #c("none", "uncertainty", "prevalence", "unc_plus_prev", "unc_plus_recs", "coverage") 
  
  # # set outpath and inputs for testing
  # dirs <- data.frame(outpath = 'broom',
  #                    inputs = 'handle')
  
  pars <- data.frame(community_file = rep(paste0(dirs$outpath, community_version, simulation_run_name, "/", community_version, sprintf(paste0("community_%i_%i_sim/", community_version, "community_%i_%i_sim_initial.rds"), n_communities, max(n_species), n_communities, max(n_species))), each = length(method)), 
                     sdm_path = rep(paste0(dirs$outpath, community_version, simulation_run_name, "/", community_version, sprintf("community_%i_%i_sim/", n_communities, max(n_species)), community_version, "species_models/"), each = length(method)), 
                     effort = as.character(paste0(dirs$inputs,"butterfly_1km_effort_layer.grd")), 
                     background = "AnnualTemp", 
                     env_data = paste0(dirs$inputs,"envdata_1km_no_corr_noNA.grd"),
                     probability_weight_adj = 1,
                     weight_adj = 1, 
                     method = method,
                     uptake = asv_version$uptake_value[asv],
                     n = 2000,
                     community_version = community_version,
                     AS_version = AS_version,
                     outPath = rep(paste0(dirs$outpath, community_version, simulation_run_name, "/", community_version, sprintf("community_%i_%i_sim/", n_communities, max(n_species))), each = length(method)))
  
  pars$rownum <- 1:nrow(pars)
  
  dim(pars)
  # head(pars)
  
  sjob <- slurm_apply(slurm_adaptive_sample, 
                      pars, 
                      jobname = paste0(community_version, '_adaptive_samp_', AS_version),
                      nodes = nrow(pars), 
                      cpus_per_node = 1, 
                      submit = TRUE,
                      slurm_options = list(partition = "short-serial-4hr",
                                           time = "01:30:00",
                                           mem = "6000",
                                           output = "sim_spp_%a.out",
                                           error = "sim_spp_%a.err",
                                           account = "short4hr"),
                      sh_template = "jasmin_submit_sh.txt")
  
}

# ### for testing
# i = 1
# 
# community_file = pars$community_file[i]
# sdm_path = pars$sdm_path[i]
# effort = pars$effort[i]
# background = pars$background[i]
# env_data =  pars$env_data[i]
# probability_weight_adj = 1
# weight_adj = 1
# method = pars$method[i]
# uptake = 0.5
# n = 2000
# community_version = pars$community_version[i]
# AS_version = pars$AS_version[i]
# outPath = pars$outPath[i]
# rownum = pars$rownum[i]
# model = c("rf", "gam", "lr")
# extent_crop = NULL