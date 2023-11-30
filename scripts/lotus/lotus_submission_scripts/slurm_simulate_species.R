
library(rslurm)

source('scripts/slurm_simulate_species_function.R')

dirs <- config::get("LOTUSpaths_sim")

# a version name that follows all the way through the community
community_version_name = 'v4'

n_communities = 12:50

pars <- data.frame(env_data = paste0(dirs$inpath, "/envdata_1km_no_corr_noNA.grd"),
                   outPath = dirs$outpath, 
                   seed = n_communities, # community number
                   max_samp = 20000, 
                   n_env = 10, 
                   n = 50,
                   det_prob = 0.2,
                   sample_across_species = TRUE,
                   effort = paste0(dirs$inpath,"butterfly_1km_effort_layer.grd"), 
                   background = "MeanDiRange",
                   community_version_name = community_version_name,
                   simulation_run_name = 'communities_1km') # the name of the run name - don't change unless changing the resolution of the area of interest.

sjob <- slurm_apply(simulate_species, pars, 
                    jobname = paste0(community_version_name, 'sim_spp'),
                    nodes = nrow(pars), 
                    cpus_per_node = 1, 
                    submit = TRUE,
                    slurm_options = list(partition = "short-serial-4hr",
                                         time = "3:59:59",
                                         mem = "10000",
                                         output = "sim_spp_%a.out",
                                         error = "sim_spp_%a.err",
                                         account = "short4hr"),
                    sh_template = "jasmin_submit_sh.txt")

