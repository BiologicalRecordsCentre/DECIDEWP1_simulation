# extract predictions and observations function
extract_predictions <- function(community_folder, community_version, AS_version, model, method, extract_preds){
  
  require(raster)
  require(readr)
  
  #read in all model files for each species
  
  model_types <- sapply(strsplit(as.character(model),","), function(x) trimws(x))
  
  models <- list.files(path = paste0(community_folder, community_version, 'species_models/'), 
                       pattern = paste0(AS_version, "_*.*(",paste(model_types, sep = "", collapse = "|"),")*.*.rdata"))
  
  # read in initial models - needed because of new naming system
  init_mods <- list.files(path = paste0(community_folder, community_version, 'species_models/'), 
                          pattern = paste0("(",paste(model_types, sep = "", collapse = "|"),")*.*initial.rdata"))
  
  # combine the two lists
  models <- c(models, init_mods)
  
  #path to initial community folder
  community_file <- paste0(community_folder, basename(as.character(community_folder)),"_initial.rds")
  
  print('these are the comm files')
  print(community_file)
  
  #read in community data
  community <- readRDS(community_file)
  
  print('loaded community file')
  
  
  #set community name
  community_name <- basename(as.character(community_folder))
  
  #average across model types for each species and method combination
  
  #get species list from length of community list
  species_list <- vector()
  for (i in 1:length(community)){species_list[i] <- paste0("Sp",i)}
  
  #get method types
  
  method_types <- as.character(sapply(strsplit(as.character(method),","), function(x) trimws(x)))
  
  eval_list <- list()
  
  # method_mod_av <- data.frame()
  # method_obsvs <- data.frame()
  
  print('starting species')
  
  #loop over species
  for (j in 1:length(species_list)){
    
    print(paste('species', j))
    
    species <- species_list[j]
    print(species)
    
    method_mod_av <- data.frame()
    method_obsvs <- data.frame()
    
    for(k in method_types){
      
      if(k == "initial"){
        models_to_read <- grep(paste0(species, "_", k, ".rdata"), models)
      } else {models_to_read <- grep(paste0(species, "_initial_AS_", k, ".rdata"), models)}
      
      if(length(models_to_read) > 0 & extract_preds){
        idx <- 1
        model_outputs <- list()
        for (l in models_to_read){
          print(models[l])
          model_output <- NULL
          try(load(paste0(community_folder, community_version, 'species_models/', models[l])))
          model_type <- model_output$model
          if(is.null(model_type)){next}#should catch corrupted rdata
          if (model_type != "rf"){
            model_preds <- model_output$predictions}
          if (model_type == "rf"){
            model_preds <- model_output$predictions[,names(model_output$predictions) %in% c("x", "y", "mean.1", "sd", "DECIDE_score.1")]
            names(model_preds) <- c("x", "y", "mean", "sd", "DECIDE_score")
          }
          model_preds <- model_preds[,1:4]#only keep mean and sd for now
          model_outputs[[idx]] <- model_preds
          names(model_outputs)[idx] <- model_type
          idx <- idx + 1
        }
        
        # average model outputs (note - not weighted by AUC currently)
        mod_average_pres_unc <- Reduce(`+`, model_outputs) / length(model_outputs)
        
        # create full dataframe of predictiosn and uncertainty
        mod_average_pres_unc <- cbind(community = community_name,
                                      method = k,
                                      species = species_list[j], 
                                      prevalence = community[[j]]$prevalence, 
                                      mod_average_pres_unc)
        
        
        # # true occurence
        # prediction <- raster::rasterFromXYZ(mod_average_pres_unc[,c('x', 'y', 'mean')])
        # true_prob_occ <- as.data.frame(raster::crop(community[[j]]$true_prob_occ, prediction), xy=T)
        # true_pa <- as.data.frame(raster::crop(community[[j]]$pres_abs, prediction), xy=T)
        
        method_mod_av <- rbind(method_mod_av, mod_average_pres_unc)
        
      } else { # this is to store information about which species and methods had 0 observations 
        
        # find and store observations from each method
        if(k == 'initial'){
          observations <- cbind(community = community_name,
                                method = k,
                                species = species_list[j], 
                                prevalence = community[[j]]$prevalence,
                                community[[j]]$observations)
        } else {
          # read in community file after AS for all AS methods
          as_comm <- readRDS(grep(AS_version, list.files(community_folder, pattern = k, full.names = T), value = T))
          observations <- cbind(community = community_name,
                                method = k,
                                species = species_list[j], 
                                prevalence = community[[j]]$prevalence,
                                as_comm[[j]]$observations)
        }
        
        method_obsvs <- rbind(method_obsvs, observations)
        
      }
      
    }#method loop
    
    dir.create(file.path(community_folder,'preds_and_obsvs'))
    
    # save model averages
    write_csv(method_mod_av, file = paste0(community_folder, 'preds_and_obsvs/', species, '_', AS_version, '_', community_name, "_model_averages.csv"))
    
    # save observations file
    write_csv(method_obsvs, file = paste0(community_folder, 'preds_and_obsvs/', species, '_', AS_version, '_', community_name, "_observations.csv"))
    
    
  }#species loop
  
  # # save model averages
  # write.csv(method_mod_av, file = paste0(community_folder, AS_version, '_', community_name, "_model_averages.csv"))
  # 
  # # save observations file
  # write.csv(method_obsvs, file = paste0(community_folder, AS_version, '_', community_name, "_observations.csv"))
  # 
  # return(list(method_mod_av, method_obsvs))
  
} #end function
