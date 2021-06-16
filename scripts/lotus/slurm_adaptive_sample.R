# function to generate new data based on existing locations and model

slurm_adaptive_sample <- function(community_file, sdm_path, model = c("rf", "gam", "lr"), method = c("none", "uncertainty", "prevalence", "unc_plus_prev", "coverage"), n = 100){
  
  #get rdata files with model outputs for each model/species (assuming communities are stored in separate folders)
  models <- list.files(path = sdm_path, pattern = paste(model, sep = "", collapse = "|"))
  
  community <- readRDS(community_file)
  
  if
  
  
  
  
}



community_file <- "N:/CEH/DECIDE/WP1_simulation/DECIDEWP1_simulation/Outputs/GB_test/_community_1000_20_sim.rds"
sdm_path <- "N:/CEH/DECIDE/WP1_simulation/DECIDEWP1_simulation/Outputs/GB_test/"
model = c("rf", "gam", "lr")