#' # Run all 10 simulated species on LOTUS
#' 
slurm_run_sim_sdm <- function(index, spdata, model, data_type, writeRas, GB, community_version, AS_version, simulation_run_name, n_communities, n_species){
  #' 
  #' ## 1. Simulate distributions (or read in simulated spp)
  library(raster)
  library(virtualspecies)
  library(dismo)
  library(tidyverse)
  library(Rfast)
  library(mgcv)
  library(randomForest)
  #library(rgdal)
  
  #load output of demo_simulatebaseline.R (could be reduced in size to speed this up)
  dirs <- config::get("LOTUSpaths")
  
  #now matches output format of slurm_simulate_species.R - can be changed?
  community <- readRDS(as.character(spdata))
  
  #' ## 2. Create input data for models
  #' 
  #' We need to extract the virtual species data we simulated and combine into a community dataset. For each species we can create a new dataset with pseudoabsences generated from the other species in the community
  #' 
  #create a pseudo-absence dataset
  source(paste0(dirs$inpath,"reformat_simulated_data.R"))
  source(paste0(dirs$inpath, "Edited_Rob_Functions.R"))
  
  #read in raster data for env data
  #read in env data frame
  
  if(GB == TRUE){
    hbv_y <- raster::stack(paste0(dirs$inpath,"envdata_1km_no_corr_noNA.grd"))
    hbv_df <- read.csv(paste0(dirs$inpath, "hbv_df_1km.csv"))} else if(GB == FALSE){
      hbv_y <- raster::stack(paste0(dirs$inpath,"hbv_y.grd")) 
      hbv_df <- readRDS(paste0(dirs$inpath, "hbv_df.rds"))
    }
  
  presences_df <- reformat_data(community, year = 2015, species_name = 'Sp')
  #head(presences_df)
  
  species_list <- unique(presences_df$species)
  
  pres_abs <- vector('list', length = length(species_list))
  
  for(s in 1:length(species_list)){
    
    pres_abs[[s]] <- cpa(spdat = presences_df, species = species_list[s], 
                         matchPres = FALSE, nAbs = 10000,
                         minYear = 2000, maxYear = 2017, recThresh = 1)
    
  }
  
  names(pres_abs) <- species_list
  
  #' Pseudoabsence weighting dependent on model, Thomas' code has this accounted for
  
  
  #' 3. Run DECIDE models 
  #' 
  #' Initially we can run just the logistic regression models as a test
  
  #source code from Thomas' workflow
  source(paste0(dirs$inpath,"getpredictions_dfsd.R"))
  
  #loop over all 10 species - set up for LOTUS
  
  #number of bootstraps
  k = 10
  
  #species index
  sp_list <- names(pres_abs)
  species <- sp_list[index]
  
  #subset envdata for species of interest  
  env_data <- subset(hbv_y, subset = community[[index]]$model_variables)
  
  #use only a 2/3 proportion of the environmental data to run the models to reduce model fit
  #env_index <- sample(1:dim(env_data_full)[3], size = round(dim(env_data_full)[3]*(2/3)), replace = FALSE)
  
  #subset environmental data for the model run
  #env_data <- env_data_full[[env_index]]
  
  #set parameters
  model <- model
  
  
  #run model for first species
  sdm <- fsdm(species = species, model = model,
              climDat = env_data, spData = pres_abs, knots_gam = 4,
              k = k, 
              write =  FALSE, outPath = paste0(dirs$outpath))
  
  #predictions
  
  preds1 <- get_predictions_dfsd(sdm, model, hbv_df)
  
  ## save files ##
  species_name <- gsub(pattern = ' ', replacement = '_', species) # get species name without space
  
  # create the output path to be in the community that the species belongs to
  outPath <- paste0(dirs$outpath, community_version, simulation_run_name,'/', community_version, "community_",n_communities,"_", n_species, "_sim/")
  
  #' Calculate very simple DECIDE score - prediction * standard deviation
  
  DECIDE_score <- preds1$mean_predictions*preds1$sd_predictions
  
  
  if(writeRas == TRUE){
    
    # create a new location to save all rasters
    dir.create(paste0(outPath, community_version, 'sdm_plots/'))
    
    # save prediction raster
    writeRaster(x = rasterFromXYZ(cbind(hbv_df$x,hbv_df$y,preds1$mean_predictions)), 
                filename = paste0(outPath, community_version, 'sdm_plots/', model, "_SDMs_", species_name, "_meanpred.grd"),
                format = 'raster', overwrite = T)
    
    # save sd raster
    writeRaster(x = rasterFromXYZ(cbind(hbv_df$x, hbv_df$y, preds1$sd_predictions)), 
                filename = paste0(outPath, community_version, 'sdm_plots/', model, "_SDMs_", species_name, "_sdpred.grd"),
                format = 'raster', overwrite = T)
    
    
    #' #' Plot maps
    #' #' 
    #' #' 
    #' png(paste0(outPath,  community_version, 'sdm_plots/', species_name,".png"), height = 200, width = 200, res = 300, units = "mm", pointsize = 14)
    #' 
    #' par(mfrow=c(3,2))
    #' par(mar = c(2,2,2,2))
    #' plot(community[[index]]$true_prob_occ, main = "Probability of occurrence")
    #' plot(community[[index]]$pres_abs, main = "Presence absence")
    #' points(community[[index]]$observations[!is.na(community[[index]]$observations$Observed),1:2], pch = 20)
    #' plot(rasterFromXYZ(cbind(hbv_df$x,hbv_df$y,preds1$mean_predictions)), main = "Predicted prob. occ")
    #' plot(rasterFromXYZ(cbind(hbv_df$x, hbv_df$y, preds1$sd_predictions)), main = "Standard deviation of predictions")
    #' plot(rasterFromXYZ(cbind(hbv_df$x, hbv_df$y, DECIDE_score)), main = "DECIDE score")
    #' 
    #' dev.off()
    
    
  }
  
  # write AUC to file for easy-access
  #write.csv(x = data.frame(raw_AUC = sdm$AUC,
  #                         meanAUC = sdm$meanAUC),
  #          file = paste0(outPath, model, "_SDMs_", species_name, "_AUC_values.csv"))
  
  # write data to file too
  #write.csv(x = sdm$Data,
  #         file = paste0(outPath, model, "_SDMs_", species_name, "_Data.csv"))
  
  # save subset model output
  # remove data from model output
  sdm$Data <- NULL
  
  community_name <- strsplit(as.character(spdata),"\\/")[[1]][10]
  
  # output of model to store
  model_output <- list(community_version,
                       AS_version,
                       community = community_name, 
                       species = species_name,
                       model = model,
                       sdm_output = lapply(sdm$Bootstrapped_models, function(x) summary(x)),
                       number_validations = k,
                       meanAUC = sdm$meanAUC,
                       predictions = data.frame(x = hbv_df$x, y = hbv_df$y, mean = preds1$mean_predictions, sd = preds1$sd_predictions, DECIDE_score = DECIDE_score))
  
  # create a new directory to store species SDMs
  dir.create(paste0(outPath, community_version, "species_models/"))
  
  print("#####     Saving output     #####") ## to check if the process is hanging on lotus
  #### Figure out how to save with the new AS version name
  save(model_output, file = paste0(outPath, community_version, "species_models/", ifelse(data_type!='initial', paste0(AS_version, '_'), ''), community_version, model, "_SDMs_GBnew_", species_name, "_", data_type, ".rdata"))
  
  print("#####     Output saved     #####")
  
  
}
