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
  env_data_full <- subset(hbv_y, subset = community[[index]]$variables)
  
  #use only a 2/3 proportion of the environmental data to run the models to reduce model fit
  env_index <- sample(1:dim(env_data_full)[3], size = round(dim(env_data_full)[3]*(2/3)), replace = FALSE)
  
  #subset environmental data for the model run
  env_data <- env_data_full[[env_index]]
  
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
    
    
    #' Plot maps
    #' 
    #' 
    png(paste0(outPath,  community_version, 'sdm_plots/', species_name,".png"), height = 200, width = 200, res = 300, units = "mm", pointsize = 14)
    
    par(mfrow=c(3,2))
    par(mar = c(2,2,2,2))
    plot(community[[index]]$true_prob_occ, main = "Probability of occurrence")
    plot(community[[index]]$pres_abs, main = "Presence absence")
    points(community[[index]]$observations[!is.na(community[[index]]$observations$Observed),1:2], pch = 20)
    plot(rasterFromXYZ(cbind(hbv_df$x,hbv_df$y,preds1$mean_predictions)), main = "Predicted prob. occ")
    plot(rasterFromXYZ(cbind(hbv_df$x, hbv_df$y, preds1$sd_predictions)), main = "Standard deviation of predictions")
    plot(rasterFromXYZ(cbind(hbv_df$x, hbv_df$y, DECIDE_score)), main = "DECIDE score")
    
    dev.off()
    
    
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

library(rslurm)

dirs <- config::get("LOTUSpaths")


## new parameters code to try and automate the parameter generation file a little more
n_species = 1:50 # vector of number of species in each community
n_communities = 1:10 # number of communities to go through
models = c('lr', 'gam', 'rf')
data_type = 'initial' # c("initial_AS_none", "initial_AS_uncertainty", "initial_AS_prevalence", "initial_AS_unc_plus_prev", "initial_AS_unc_plus_recs", "initial_AS_coverage") # 'initial'

# name of the versions we are running - so we're not overwriting things
# one for community-level which includes the community folders and species models folders
community_version = 'v3'

# One name for the adaptive sampling round to allow us to create different sampling methods of the same initial community
# This doesn't get used if running only the initial model, but does get used when running the adaptive sampling methods.
AS_version = 'asv1'

# the name of the simulation run - same as slurm_simulate species
simulation_run_name = 'communities_1km'

# # set commpath to something for testing (then delete 'dirs$')
# dirs <- data.frame(commpath = 'blob')

# pars data frame
pars <- data.frame(index = rep(n_species, length(n_communities)*length(models)*length(data_type)),
                   spdata = rep(sprintf(
                     paste0(dirs$commpath, community_version, simulation_run_name,"/", community_version,
                            "community_%i_%i_sim/", ifelse(data_type!='initial', paste0(AS_version, '_'), ''), community_version, "community_%i_%i_sim_%s.rds"), 
                     rep(rep(n_communities, each = length(models)), each = length(data_type)), max(n_species), 
                     rep(rep(n_communities, each = length(models)), each = length(data_type)), max(n_species), data_type
                   ), each = length(n_species)),
                   model = rep(rep(rep(models, length(n_communities)), each = length(data_type)), each = length(n_species)),
                   data_type = rep(rep(data_type, length(n_communities)*length(models)), each = length(n_species)), 
                   writeRas = FALSE,
                   GB = TRUE,
                   community_version = community_version,
                   simulation_run_name = simulation_run_name,
                   AS_version = AS_version,
                   n_communities = rep(rep(rep(n_communities, each = length(models)), each = length(data_type)), each = length(n_species)),
                   n_species = max(n_species)
)

dim(pars)

#test with subset of runs
#pars <- pars[-c(1,78,103,166,208,294,315,352,422, 484,517, 569, 646,675, 735),]

# # resubmit long runs - job numbers from lotus
# resub_rows <- c(1486, 1485, 1455, 1456, 1388, 1389, 1395, 1396, 1365, 1366, 1335, 1336, 1305, 1306, 1245, 1246, 1238, 1239)+1
# pars <- pars[resub_rows,]

#### slurm apply call
sdm_slurm <- slurm_apply(slurm_run_sim_sdm,
                         params = pars,
                         jobname = paste0(community_version, '_sdm_simulated_species'),
                         nodes = length(pars$index),
                         cpus_per_node = 1,
                         slurm_options = list(partition = 'short-serial',#-4hr',
                                              time = '23:59:59',
                                              mem = 10000,
                                              output = "sim_sdm_%a.out",
                                              error = "sim_sdm_%a.err"),#,
                         # account = "short4hr"),
                         sh_template = "jasmin_submit_sh.txt",
                         submit = T)
pars$BatchID <- sdm_slurm$jobid
pars$JobID <- 0:(nrow(pars)-1)#slurm job ID
write.csv(pars, paste0('_rslurm_', community_version, '_sdm_simulated_species/pars.csv'))#to match to error files


