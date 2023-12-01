## percentage increase in data

### run on lotus!!

dirs <- config::get("LOTUSpaths")

# name of the community version we are running - so we're not overwriting things, keep same as for slurm_run_sim_sdm
community_version = 'v4'

# name of the adaptive sampling version we are looking to evaluate
AS_version = 1:4

# the name of the simulation run - same as slurm_simulate species
simulation_run_name = 'communities_1km'

n_communities = 1:50

n_species = 1:50

community_folder = paste0(dirs$outpath, community_version, simulation_run_name, "/", community_version, sprintf("community_%i_%i_sim/", n_communities, max(n_species)))

method = c('none', 'uncertainty', 'prevalence', 'unc_plus_prev', 'unc_plus_recs', 'coverage')

observations_out <- data.frame()

# community loop
for(com in n_communities) {
  print(com)
  
  inital_comm <- readRDS(list.files(community_folder[com], pattern = 'initial.rds', full.names = TRUE))
  
  # AS loop
  for(asv in AS_version) {
    
    asv_oi <- list.files(community_folder[com], pattern = paste0('asv', asv), full.names = TRUE)
    
    for(meth in 1:length(method)) {
      
      as_file <- readRDS(grep(method[meth], asv_oi, value = T))
      
      for(spec in n_species) {
        
        init_obvs <- sum(inital_comm[[spec]]$observations$Observed, na.rm = T)    
        new_observations <- sum(as_file[[spec]]$observations$Observed, na.rm = T)
        
        observations_out <- rbind(observations_out, 
                                  data.frame(community_version = community_version, 
                                             community = com,
                                             as_version = AS_version,
                                             species = spec,
                                             meth = method[meth],
                                             initial_obvs = init_obvs,
                                             new_observations = new_observations - init_obvs,
                                             perc_increase = (new_observations - init_obvs)/init_obvs*100))
        
      }
      
    }
    
  }
  
}


write.csv(observations_out, file = paste0(dirs$outpath, community_version, simulation_run_name, 
                             "/", community_version, "n_new_obvs_perc_increase_", min(n_communities), "_", max(n_communities), "_spp", max(n_species), ".csv"))
