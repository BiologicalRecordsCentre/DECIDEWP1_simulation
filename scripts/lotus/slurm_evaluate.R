slurm_evaluate <- function(community_folder, model, method){
  
  #read in all model files for each species
  
  model_types <- sapply(strsplit(model,","), function(x) trimws(x))
  
  models <- list.files(path = community_folder, pattern = paste0("(",paste(model_types, sep = "", collapse = "|"),")*.*.rdata"))
  
  #identify community file based on naming structure
  community_file <- paste0(community_folder, basename(community_folder),".rds")
  
  #read in community data
  community <- readRDS(community_file)
  
  #average across model types for each species and method combination
  
  #get species list from length of community list
  species_list <- vector()
  for (i in 1:length(community)){species_list[i] <- paste0("Sp",i)}
  
  #get method types
  
  method_types <- as.character(sapply(strsplit(method,","), function(x) trimws(x)))
  
  eval_list <- list()
  
  #loop over species
  for (j in 1:length(species_list)){
    
    species <- species_list[j]
    
    method_eval <- data.frame()
    for(k in method_types){
    
      if(k == "initial"){
        models_to_read <- grep(paste0(species, "_", k), models)
      } else {models_to_read <- grep(paste0(species, "_AS_*(", k, ")"), models)}
    
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
    auc <- as.numeric(pROC::auc(eval_df$true_pa, eval_df$prediction))
    
    method_eval <- rbind(method_eval, data.frame(method = k, mse = mse, corr = corr, auc = auc, species = species))
  
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
      prevalence[j] <- sum(getValues(community[[j]]$pres_abs), na.rm=TRUE)/nrow(mod_average)
    }
  } else {prevalence <- sapply(community, function(x) x$prevalence)}
  
  eval_table$prevalence <- prevalence[as.numeric(sapply(strsplit(eval_table$species, split = "Sp"), function(x) x[[2]]))]
  
  return(eval_table)
  
} #end function


library(rslurm)

dirs <- config::get("LOTUSpaths")

## index file
pars <- data.frame(community_folder = paste0(dirs$commpath, "community_1_50_sim/"), model = "rf, gam, lr", method = "initial, none, uncertainty, prevalence, unc_plus_prev, coverage")

#### slurm apply call
sdm_slurm <- slurm_apply(slurm_evaluate,
                         params = pars,
                         jobname = 'evaluate',
                         nodes = length(pars$community_folder),
                         cpus_per_node = 1,
                         slurm_options = list(partition = 'test',
                                              time = '0:14:59',
                                              mem = 3000,
                                              output = "sim_eval_%a.out",
                                              error = "sim_eval_%a.err"),
                         sh_template = "jasmin_submit_sh.txt",
                         submit = T)

