slurm_evaluate <- function(community_folder, community_version, AS_version, model, method){
  
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
  
  #read in community data
  community <- readRDS(community_file)
  
  #set community name
  community_name <- basename(as.character(community_folder))
  
  #average across model types for each species and method combination
  
  #get species list from length of community list
  species_list <- vector()
  for (i in 1:length(community)){species_list[i] <- paste0("Sp",i)}
  
  #get method types
  
  method_types <- as.character(sapply(strsplit(as.character(method),","), function(x) trimws(x)))
  
  eval_list <- list()
  
  #loop over species
  for (j in 1:length(species_list)){
    
    species <- species_list[j]
    
    method_eval <- data.frame()
    for(k in method_types){
      
      if(k == "initial"){
        models_to_read <- grep(paste0(species, "_", k, ".rdata"), models)
      } else {models_to_read <- grep(paste0(species, "_initial_AS_", k, ".rdata"), models)}
      
      if(length(models_to_read) > 0){
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
        
        #average model outputs (note - not weighted by AUC currently)
        mod_average <- Reduce(`+`, model_outputs) / length(model_outputs)
        
        #extract basic metrics
        
        
        prediction <- raster::rasterFromXYZ(mod_average[,1:3])
        true_prob_occ <- raster::crop(community[[j]]$true_prob_occ, prediction)
        true_pa <- raster::crop(community[[j]]$pres_abs, prediction)
        
        mse <- mean((raster::getValues(true_prob_occ)-raster::getValues(prediction))^2, na.rm=TRUE)
        medianse <- median((raster::getValues(true_prob_occ)-raster::getValues(prediction))^2, na.rm=TRUE)
        corr <- cor(raster::getValues(true_prob_occ), raster::getValues(prediction), use = "pairwise.complete")
        auc <- as.numeric(pROC::auc(raster::getValues(true_pa), raster::getValues(prediction), quiet = TRUE))
        
        # extract mean and sd error
        mabse <- mean((raster::getValues(true_prob_occ)-raster::getValues(prediction)), na.rm=TRUE)
        sdabse <- sd((raster::getValues(true_prob_occ)-raster::getValues(prediction)), na.rm=TRUE)
        
        occprobe <- sum(raster::getValues(true_prob_occ), na.rm = T) - sum(raster::getValues(prediction), na.rm=TRUE)
        
        #extract mean and max stdev
        mean_sd <- mean(mod_average$sd, na.rm = TRUE)
        max_sd <- max(mod_average$sd, na.rm = TRUE)
        
        method_eval <- rbind(method_eval, data.frame(method = k, mse = mse, medianse = medianse, corr = corr, auc = auc, mean_sd = mean_sd, max_sd = max_sd, species = species, mabse = mabse, sdabse = sdabse, occprobe = occprobe))
        
      }
    }#method loop
    
    if(nrow(method_eval) >0 ){
      eval_list[[j]] <- method_eval} else {eval_list[[j]] <- NULL}
    
  }#species loop
  
  
  eval_table <- do.call("rbind", eval_list)
  
  #extract prevalence values for all species - calculate if not in community object
  if(is.null(community[[1]]$prevalence)){
    prevalence <- vector()
    for (j in 1:length(species_list)){
      prevalence[j] <- sum(raster::getValues(community[[j]]$pres_abs), na.rm=TRUE)/nrow(mod_average)
    }
  } else {prevalence <- sapply(community, function(x) x$prevalence)}
  
  eval_table$prevalence <- prevalence[as.numeric(sapply(strsplit(as.character(eval_table$species), split = "Sp"), function(x) x[[2]]))]
  
  eval_table$community <- community_name
  
  write.csv(eval_table, file = paste0(community_folder, community_name, "_evaluation_table2.csv"))
  
  ### different format
  init_tab <- eval_table[eval_table$method =='initial',]
  colnames(init_tab) <- paste0('initial_', colnames(init_tab))
  et <- eval_table
  
  et$init_mse <- init_tab$initial_mse[match(paste0(et$species, et$community),
                                            paste0(init_tab$initial_species, init_tab$initial_community))]
  et$init_corr <- init_tab$initial_corr[match(paste0(et$species, et$community),
                                              paste0(init_tab$initial_species, init_tab$initial_community))]
  et$init_auc <- init_tab$initial_auc[match(paste0(et$species, et$community),
                                            paste0(init_tab$initial_species, init_tab$initial_community))]
  
  # alternate format
  write.csv(et, file = paste0(community_folder, AS_version, '_', community_name, "_evaluation_table_alt2.csv"))
  
  
} #end function
