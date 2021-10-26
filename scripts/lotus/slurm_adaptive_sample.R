## script to submit adaptive sampling

# source function to run
source("scripts/slurm_adaptive_sample_function.R")

library(rslurm)

dirs <- config::get("LOTUSpaths_AS")

# name of the version we are running - so we're not overwriting things, keep same as for slurm_run_sim_sdm
version_name = 'v1'

# the name of the simulation run - same as slurm_simulate species
simulation_run_name = 'communities_1km'

# number ofc communities
n_communities = 1

# number of species in each community - used only in the parameter file to allow runs with different numbers of species
n_species = 50

method = c("none", "uncertainty", "prevalence", "unc_plus_prev", "unc_plus_recs", "coverage")

pars <- data.frame(community_file = rep(paste0(dirs$outpath, version_name, simulation_run_name, "/", version_name, sprintf(paste0("community_%i_%i_sim/", version_name, "community_%i_%i_sim_initial.rds"), n_communities, max(n_species), n_communities, max(n_species))), each = length(method)), 
                   sdm_path = rep(paste0(dirs$outpath, version_name, simulation_run_name, "/", version_name, sprintf("community_%i_%i_sim/", n_communities, max(n_species)), version_name, "species_models/"), each = length(method)), 
                   effort = paste0(dirs$inputs,"butterfly_1km_effort_layer.grd"), 
                   background = "AnnualTemp", 
                   env_data = paste0(dirs$inputs,"envdata_1km_no_corr_noNA.grd"), 
                   weight_adj = 1, 
                   method = method, 
                   n = 2000,
                   version_name = version_name,
                   outPath = rep(paste0(dirs$outpath, version_name, simulation_run_name, "/", version_name, sprintf("community_%i_%i_sim/", n_communities, max(n_species))), each = length(method)))


sjob <- slurm_apply(slurm_adaptive_sample, pars, 
                    jobname = 'adaptive_samp_test',
                    nodes = nrow(pars), 
                    cpus_per_node = 1, 
                    submit = TRUE,
                    slurm_options = list(partition = "test", # "short-serial-4hr",
                                         time = "00:03:59",
                                         mem = "6000",
                                         output = "sim_spp_%a.out",
                                         error = "sim_spp_%a.err"),
                                         # account = "short4hr"),
                    sh_template = "jasmin_submit_sh.txt")

