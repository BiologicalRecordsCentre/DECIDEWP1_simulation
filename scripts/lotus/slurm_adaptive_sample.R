# function to generate new data based on existing locations and model

slurm_adaptive_sample <- function(community_file, sdm_path, model = c("rf", "gam", "lr"), method = c("none", "uncertainty", "prevalence", "unc_plus_prev", "coverage"), n = 100){
  
  #get rdata files with model outputs for each model/species (assuming communities are stored in separate folders)
  models <- list.files(path = sdm_path, pattern = paste(model, sep = "", collapse = "|"))
  
  community <- readRDS(community_file)
  
  #get species list from length of community list
  species_list <- vector()
  for (i in 1:length(community)){species_list[i] <- paste0("Sp",i)}
  
  #for each species on the list, extract the relevant model outputs (this allows for some models to fail for some species)
  
  community_preds <- list()
  
  for (j in 1:length(species_list)){
  
    species <- species_list[j]
    #get model outputs for this species
    models_to_read <- grep(species, models)
  
    #load outputs into list and extract prediction table
    model_outputs <- list()
    idx <- 1
    for (k in models_to_read){
      load(paste0(sdm_path, models[k]))
      model_type <- model_output$model
      if (model_type != "rf"){
        model_preds <- model_output$predictions}
      if (model_type == "rf"){
        model_preds <- model_output$predictions[,names(model_output$predictions) %in% c("x", "y", "mean.1", "sd", "DECIDE_score.1")]
        names(model_preds) <- c("x", "y", "mean", "sd", "DECIDE_score")
      }
      model_outputs[[idx]] <- model_preds
      names(model_outputs)[idx] <- model_type
      idx <- idx + 1
    }
    
    #average model outputs (note - not weighted by AUC currently)
    mod_average <- Reduce(`+`, model_outputs) / length(model_outputs)
    
    #store only the model average for now - could edit to store the individual model outputs if needed
    community_preds[[j]] <- mod_average
    names(community_preds)[j] <- species
  
  }
  
  #average across all species in community to obtain a single prevalence, uncertainty and DECIDE score 
  
  community_scores <- Reduce(`+`, community_preds)/length(community_preds)
  
  if (method == "none"){}
  
  if (method == "uncertainty"){}
  
  if (method == "prevalence"){}
  
  if (method == "unc_plus_prev"){
    cell_weights <- community_scores$DECIDE_score/sum(community_scores$DECIDE_score, na.rm=TRUE)
    #assign NA values the average weight
    cell_weights[is.na(cell_weights)] <- mean(cell_weights, na.rm= TRUE)
    #sample new locations according to cell weights
    new_locs <- sample(1:nrow(community_scores), size = n, replace = FALSE, prob = cell_weights)
    new_coords <- community_scores[new_locs, 1:2]
  }
  
  if (method == "coverage"){}
  
}



community_file <- "N:/CEH/DECIDE/WP1_simulation/DECIDEWP1_simulation/Outputs/GB_test/_community_1000_20_sim.rds"
sdm_path <- "N:/CEH/DECIDE/WP1_simulation/DECIDEWP1_simulation/Outputs/GB_test/"
model = c("rf", "gam", "lr")