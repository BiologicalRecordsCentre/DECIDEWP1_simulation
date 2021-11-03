
## Combine communities from lotus

library(rslurm)

dirs <- config::get("LOTUSpaths")

# name of the version we are running - so we're not overwriting things, keep same as for slurm_run_sim_sdm
version_name = 'v2'

# the name of the simulation run - same as slurm_simulate species
simulation_run_name = 'communities_1km'

n_communities = 1:10

n_species = 1:50

## index file
community_folder = paste0(dirs$outpath, version_name, simulation_run_name, "/", version_name, sprintf("community_%i_%i_sim/", n_communities, max(n_species)))

out_files <- list()

for(i in 1:length(community_folder)) {
  
  f_to_read <- list.files(path = community_folder[i], pattern = "_evaluation_table_alt.csv", full.names = TRUE)
  
  if(length(f_to_read)==0) next
  
  out_files[[i]] <- read.csv(f_to_read)
  
}

out <- do.call(rbind, out_files)

write.csv(out, file = paste0(dirs$outpath, version_name, simulation_run_name, 
                             "/", version_name, "combined_outputs_comm", min(n_communities), "_", max(n_communities), "_spp", max(n_species), ".csv"))
          
