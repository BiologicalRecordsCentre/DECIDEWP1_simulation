
library(rslurm)

dirs <- config::get("LOTUSpaths_AS")

# name of the versions we are running - so we're not overwriting things
# Three different versions, one for community-level which includes the community folders and species models folders
community_version = 'v2'

# # one for the initial model - will almost certainly always be the same as the community_folder - currently not included
# # because the community is essentially only the initial .rds file 
# initial_community_version = 'v2'

# and an adaptive sampling version, which is if we want to run the adaptive sampling 
# process more than once - these outputs are stored in the same place as the old outputs
# must always be prefixed by asv
AS_version = 'asv1'

# the name of the simulation run - same as slurm_simulate species
simulation_run_name = 'communities_1km'

# number ofc communities
n_communities = 1:10

# number of species in each community - used only in the parameter file to allow runs with different numbers of species
n_species = 1:50

# the adaptive sampling methods to use 
method = c("none", "uncertainty", "prevalence", "unc_plus_prev", "unc_plus_recs", "coverage") 

# community files
community_file = (paste0(dirs$outpath, community_version, simulation_run_name, "/", community_version, sprintf(paste0("community_%i_%i_sim/", community_version, "community_%i_%i_sim_initial.rds"), n_communities, max(n_species), n_communities, max(n_species))))

out_df <- list()

# loop through community files
for(coms in 1:length(community_file)) {
  
  print(coms)
  
  comm_file <- readRDS(community_file[coms])
  
  sp_df <- data.frame(community_name = rep(NA, 50),
                       species_name = rep(NA, 50),
                       real_locs_nums = rep(NA, 50),
                       observed_locs_nums = rep(NA, 50))
  
  # loop through species files
  for(spec in 1:length(comm_file)) {
    
    sp_df$community_name[spec] <- paste0('Community_', coms)
    sp_df$species_name[spec] <- paste0('Sp',spec)
    sp_df$real_locs_nums[spec] <- sum(comm_file[[spec]]$observations$Real)
    sp_df$observed_locs_nums[spec] <- sum(comm_file[[spec]]$observations$Observed, na.rm=T)
    
  }
  
  out_df[[coms]] <- sp_df
  
}

comm_df <- do.call('rbind', out_df)

write.csv(comm_df, file = paste0(dirs$outpath, community_version, simulation_run_name, 
                                 "/", AS_version, "_", community_version, "_number_of_observations_", min(n_communities), "_", max(n_communities), "_spp", max(n_species), ".csv"))


