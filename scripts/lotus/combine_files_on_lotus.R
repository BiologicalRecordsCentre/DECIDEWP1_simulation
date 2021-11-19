
## Combine communities from lotus

library(rslurm)

dirs <- config::get("LOTUSpaths")

# the name of the simulation run - same as slurm_simulate species
simulation_run_name = 'communities_1km'

# name of the community version we are running - so we're not overwriting things, keep same as for slurm_run_sim_sdm
community_version = 'v2'

# name of the adaptive sampling version we are looking to evaluate
AS_version = 'asv1'

n_communities = 1:10

n_species = 1:50

## index file
community_folder = paste0(dirs$outpath, community_version, simulation_run_name, "/", community_version, sprintf("community_%i_%i_sim/", n_communities, max(n_species)))

out_files <- list()

for(i in 1:length(community_folder)) {
  
  f_to_read <- list.files(path = community_folder[i], pattern = paste0(AS_version, "*.*_evaluation_table_alt.csv"), full.names = TRUE)
  
  if(length(f_to_read)==0) next
  
  out_files[[i]] <- read.csv(f_to_read)
  
}

out <- do.call(rbind, out_files)

write.csv(out, file = paste0(dirs$outpath, community_version, simulation_run_name, 
                             "/", AS_version, "_", community_version, "combined_outputs_comm", min(n_communities), "_", max(n_communities), "_spp", max(n_species), ".csv"))
          
