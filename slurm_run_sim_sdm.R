#' # Run all 10 simulated species on LOTUS
#' 
slurm_run_sim_sdm <- function(index){
#' 
  #' ## 1. Simulate distributions (or read in simulated spp)
  library(raster)
  library(virtualspecies)
  library(dismo)
  library(rgdal)
  
  #load output of demo_simulatebaseline.R (could be reduced in size to speed this up)
  load("virt_comm_10spp.Rdata")
  
  #' ## 2. Create input data for models
  #' 
  #' We need to extract the virtual species data we simulated and combine into a community dataset. For each species we can create a new dataset with pseudoabsences generated from the other species in the community
  #' 
  #create a pseudo-absence dataset
  pa_sets <- list()
  
  for(i in 1:length(sp.obs)){
    #extract occurrence records
    occs <- sp.obs[[i]]$sample.points[!is.na(sp.obs[[i]]$sample.points$Observed),]
    abs <- data.frame()
    for (j in 1:length(sp.obs)){
      if(j != i) {
        #observations of other species treated as absences for species of interest
        abs <- rbind(abs, sp.obs[[j]]$sample.points[!is.na(sp.obs[[j]]$sample.points$Observed),])
      }
    }
    abs$Observed <- 0
    #check no duplicates
    pa_data <- rbind(occs, abs)
    dups <- pa_data[duplicated(pa_data[1:2]),]
    #always remove duplicates - because occurrences bound first then any absences in same location as occurrences will be removed, also removed any duplicates in pseudo absence data
    if(nrow(dups) == 0) {next} else {
      pa_data <- pa_data[!duplicated(pa_data[1:2]),]
    }
    pa_sets[[i]] <- list(pa_data = pa_data, Presence = pa_data[pa_data$Observed == 1,1:2], pseudoAbsence = pa_data[pa_data$Observed == 0,1:2])
    names(pa_sets)[i] <- paste0("Sp",i)
  }
  
  #' Pseudoabsence weighting dependent on model, Thomas' code has this accounted for
  
  
  #' 3. Run DECIDE models 
  #' 
  #' Initially we can run just the logistic regression models as a test
  
  #source Edited Rob Functions
  source("Edited_Rob_Functions.R")
  
  #source code from Thomas' workflow
  source("getpredictions.R")
  
  #read in raster data for env data
  hbv_y <- raster::stack("hbv_y.grd") 
  
  #loop over all 10 species - set up for LOTUS
  
  #number of bootstraps
  k = 10
  
  #species index
  sp_list <- names(pa_sets)
  species <- sp_list[index]
    
  #subset envdata for species of interest  
  env_data <- subset(hbv_y, subset = virt_comm1[[index]]$details$variables)
  
  #set parameters
  model <- "lr"
  
    
  #run model for first species
  sdm <- fsdm(species = species, model = model,
                climDat = env_data, spData = pa_sets, knots_gam = -1,
                k = k, 
                write =  TRUE, outPath = "lr_outs/")
    
  #predictions
    
  preds1 <- get_predictions(sdm, model, env_data)
    
  ## save files ##
  species_name <- gsub(pattern = ' ', replacement = '_', species) # get species name without space
    
  outPath <- paste0(getwd(), "/Outputs/")
    
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
    
  DECIDE_score <- preds1$mean_predictions*preds1$quant_range
   
  writeRaster(x = DECIDE_score, 
                  filename = paste0(outPath, model, "_SDMs_", species_name, "_DECIDEscore.grd"),
                  format = 'raster', overwrite = T)
      
  #' Plot maps
  #' 
  png(paste0("Plots/", species_name,".png"), height = 200, width = 200, res = 300, units = "mm", pointsize = 14)
      
  par(mfrow=c(3,2))
  par(mar = c(2,2,2,2))
  plot(virt_comm1[[index]], main = "Environmental suitability")
  plot(pa[[index]]$probability.of.occurrence, main = "Probability of occurrence")
  plot(pa[[index]]$pa.raster, main = "Presence absence")
  points(sp.obs[[index]]$sample.points[is.na(sp.obs[[index]]$sample.points$Observed),1:2], pch = 20)
  plot(preds1$mean_predictions, main = "Predicted prob. occ")
  plot(preds1$quant_range, main = "Quantile range of predictions")
  plot(DECIDE_score, main = "DECIDE score")
      
  dev.off()
  }
  
