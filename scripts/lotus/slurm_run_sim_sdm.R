#' # Run all 10 simulated species on LOTUS
#' 
slurm_run_sim_sdm <- function(index, spdata){
  #' 
  #' ## 1. Simulate distributions (or read in simulated spp)
  library(raster)
  library(virtualspecies)
  library(dismo)
  library(tidyverse)
  #library(rgdal)
  
  #load output of demo_simulatebaseline.R (could be reduced in size to speed this up)
  dirs <- config::get("LOTUSpaths")
  
  #now matches output format of slurm_simulate_species.R - can be changed?
  community <- readRDS(spdata)[[1]]
  
  #' ## 2. Create input data for models
  #' 
  #' We need to extract the virtual species data we simulated and combine into a community dataset. For each species we can create a new dataset with pseudoabsences generated from the other species in the community
  #' 
  #create a pseudo-absence dataset
  source(paste0(dirs$inpath,"reformat_simulated_data.R"))
  source(paste0(dirs$inpath, "Edited_Rob_Functions.R"))
  
  presences_df <- reformat_data(community, year = 2015, species_name = 'Sp')
  head(presences_df)
  
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
  source(paste0(dirs$inpath,"getpredictions.R"))
  
  #read in raster data for env data
  hbv_y <- raster::stack(paste0(dirs$inpath,"hbv_y.grd")) 
  
  #loop over all 10 species - set up for LOTUS
  
  #number of bootstraps
  k = 10
  
  #species index
  sp_list <- names(pres_abs)
  species <- sp_list[index]
  
  #subset envdata for species of interest  
  env_data <- subset(hbv_y, subset = community[[index]]$variables)
  
  #set parameters
  model <- "lr"
  
  
  #run model for first species
  sdm <- fsdm(species = species, model = model,
              climDat = env_data, spData = pres_abs, knots_gam = -1,
              k = k, 
              write =  TRUE, outPath = paste0(dirs$outpath,"lr_outs/"))
  
  #predictions
  
  preds1 <- get_predictions(sdm, model, env_data)
  
  ## save files ##
  species_name <- gsub(pattern = ' ', replacement = '_', species) # get species name without space
  
  outPath <- dirs$outpath
  
  # save prediction raster
  writeRaster(x = preds1$mean_predictions, 
              filename = paste0(outPath, model, "_SDMs_", species_name, "_meanpred.grd"),
              format = 'raster', overwrite = T)
  
  # save quantile max min
  writeRaster(x = preds1$quant_minmax, 
              filename = paste0(outPath, model, "_SDMs_", species_name, "_quantilemaxmin.grd"),
              format = 'raster', overwrite = T)
  
  # save quantile range raster
  writeRaster(x = preds1$quant_range, 
              filename = paste0(outPath, model, "_SDMs_", species_name, "_quantilerange.grd"),
              format = 'raster', overwrite = T)
  
  # write AUC to file for easy-access
  write.csv(x = data.frame(raw_AUC = sdm$AUC,
                           meanAUC = sdm$meanAUC),
            file = paste0(outPath, model, "_SDMs_", species_name, "_AUC_values.csv"))
  
  # write data to file too
  write.csv(x = sdm$Data,
            file = paste0(outPath, model, "_SDMs_", species_name, "_Data.csv"))
  
  # save subset model output
  # remove data from model output
  sdm$Data <- NULL
  
  # output of model to store
  model_output <- list(species = species_name,
                       model = model,
                       sdm_output = sdm,
                       number_validations = k)
  
  save(model_output, file = paste0(outPath, model, "_SDMs_", species_name, 
                                   ".rdata"))
  
  
  #' Calculate very simple DECIDE score - prediction * quantile range
  
  DECIDE_score <- preds1$mean_predictions*sqrt(preds1$quant_range)
  
  writeRaster(x = DECIDE_score, 
              filename = paste0(outPath, model, "_SDMs_", species_name, "_DECIDEscore.grd"),
              format = 'raster', overwrite = T)
  
  #' Plot maps
  #' 
  png(paste0(dirs$outpath,"Plots/", species_name,".png"), height = 200, width = 200, res = 300, units = "mm", pointsize = 14)
  
  par(mfrow=c(3,2))
  par(mar = c(2,2,2,2))
  plot(community[[index]]$true_prob_occ, main = "Probability of occurrence")
  plot(community[[index]]$pres_abs, main = "Presence absence")
  points(community[[index]]$observations[!is.na(community[[index]]$observations$Observed),1:2], pch = 20)
  plot(preds1$mean_predictions, main = "Predicted prob. occ")
  plot(preds1$quant_range, main = "Quantile range of predictions")
  plot(DECIDE_score, main = "DECIDE score")
  
  dev.off()
}

## index file
pars <- data.frame(index = seq(1:20), spdata = "/gws/nopw/j04/ceh_generic/susjar/DECIDE/_rslurm_sim_spp/results_0.RDS") # number of species

library(rslurm)

dirs <- config::get("LOTUSpaths")

#### slurm apply call
sdm_slurm <- slurm_apply(slurm_run_sim_sdm,
                         params = pars,
                         jobname = 'sdm_simulated_species',
                         nodes = length(pars$index),
                         cpus_per_node = 1,
                         slurm_options = list(partition = 'short-serial',
                                              time = '23:59:59',
                                              mem = 20000,
                                              output = "sim_sdm_%a.out",
                                              error = "sim_sdm_%a.err"),
                         sh_template = "jasmin_submit_sh.txt",
                         submit = T)


