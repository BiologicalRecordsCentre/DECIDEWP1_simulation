slurm_evaluate <- function(community_folder, model, method){
  
  #read in all model files except initial for each species
  
  models <- list.files(path = community_folder, pattern = paste0("(",paste(model, sep = "", collapse = "|"),")*AS.*.rdata"))
  
  #identify community file based on naming structure
  community_file <- paste0(community_folder, basename(community_folder),".rds")
  
  #read in community data
  community <- readRDS(community_file)
  
  #average across model types for each species and method combination
  
  #get species list from length of community list
  species_list <- vector()
  for (i in 1:length(community)){species_list[i] <- paste0("Sp",i)}
  
  #loop over species
  for (j in 1:length(species_list)){
    
    species <- species_list[j]
    
    method_eval <- data.frame()
    for(k in method){
    
    models_to_read <- grep(paste0(species, "_AS_*(", k, ")"), models)
    
    if(length(models_to_read > 0)){
      idx <- 1
      model_outputs <- list()
      for (l in models_to_read){
        model_output <- NULL
        try(load(paste0(community_folder, models[l])))
        model_type <- model_output$model
        if (model_type != "rf"){
          model_preds <- model_output$predictions}
        if (model_type == "rf"){
          model_preds <- model_output$predictions[,names(model_output$predictions) %in% c("x", "y", "mean.1", "sd", "DECIDE_score.1")]
          names(model_preds) <- c("x", "y", "mean", "sd", "DECIDE_score")
                    }
        model_preds <- model_preds[,1:3]#only keep mean for now
        model_outputs[[idx]] <- model_preds
        names(model_outputs)[idx] <- model_type
        idx <- idx + 1
        }
    
    #average model outputs (note - not weighted by AUC currently)
    mod_average <- Reduce(`+`, model_outputs) / length(model_outputs)
      
    #extract basic metrics
    
    true_prob_occ <- raster::as.data.frame(community[[j]]$true_prob_occ, xy=TRUE)
    true_pa <- raster::as.data.frame(community[[j]]$pres_abs, xy=TRUE)
    prediction <- mod_average
    
    eval_df <- Reduce(function(df1, df2) merge(df1, df2, by = c("x", "y")), list(true_prob_occ, true_pa, prediction))
    names(eval_df) <- c("x", "y", "true_prob_occ", "true_pa", "prediction")
    
    mse <- mean((eval_df$true_prob_occ-eval_df$prediction)^2, na.rm=TRUE)
    corr <- cor(eval_df$true_prob_occ, eval_df$prediction)
    
    method_eval <- rbind(method_eval, data.frame(method = k, mse = mse, corr = corr))
      }
    }#method loop
    
    
  }#species loop
  
  
  
  
  
} #end function


community_folder <- "N:/CEH/DECIDE/WP1_simulation/DECIDEWP1_simulation/Outputs/community_1_50_sim/"
model = c("rf", "gam", "lr")
method = c("none", "uncertainty", "prevalence", "unc_plus_prev", "coverage")
