# function to generate new data based on existing locations and model

slurm_adaptive_sample <- function(rownum, community_file, sdm_path, effort, background, env_data, extent_crop = NULL, probability_weight_adj, weight_adj, model = c("rf", "gam", "lr"), method, n = 100, uptake = NULL, community_version, AS_version, outPath){
  
  # print the row number from the pars file
  print(rownum)
  print(method)
  print(community_file)
  
  #get rdata files with model outputs for each model/species (assuming communities are stored in separate folders) - only read initial models
  models <- list.files(path = as.character(sdm_path), pattern = paste0("(",paste(model, sep = "", collapse = "|"),")*initial.rdata"))
  
  # import simulated community data
  community <- readRDS(as.character(community_file))
  
  #extract prevalence vector
  prevalence_vec <- sapply(community, function(x) x$prevalence)
  
  #import env_data if specified
  if(!is.null(env_data)){
    env <- raster::stack(as.character(env_data))
    
    #crop to extent if specified
    if(!is.null(extent_crop)){
      e <- as(extent_crop, "SpatialPolygons")
      sp::proj4string(e) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
      
      e.geo <- sp::spTransform(e, CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs"))
      
      env_extent <- raster::crop(env, e.geo)
      
    } else {env_extent <- env}
    
  }
  
  #set background if given, can indicate a layer in env_data or be a filepath to a raster
  if(is.numeric(background)){
    bg_layer <- env_extent[[background]]
  } else if ((is.character(background)|is.factor(background)) & !grepl("\\.", background)) {bg_layer <- raster::subset(env_extent, background)} else if ((is.character(background)|is.factor(background)) & grepl("\\.", background)) {bg_layer <- raster::raster(as.character(background))} else {bg_layer <- NULL}
  
  #extract effort layer from raster if provided (note currently uses layers in existing raster stack, could read in other layers)
  if(is.numeric(effort)){eff_layer <- env_extent[[effort]]} else if((is.character(effort)|is.factor(effort)) & !grepl("\\.", effort)) {eff_layer <- raster::subset(env_extent,effort)} else if ((is.character(effort)|is.factor(effort)) & grepl("\\.", effort)) {eff_layer <- raster::raster(as.character(effort))} else  {eff_layer <- NULL}
  
  if(is.null(eff_layer)){eff_weights <- (env_extent[[1]]*0)+1} else if (is.null(bg_layer)){
    eff_weights <- eff_layer/weight_adj} else {eff_weights <- (bg_layer/bg_layer) + (eff_layer/weight_adj)}
  
  eff_df <- raster::as.data.frame(eff_weights, xy=TRUE, na.rm=TRUE) 
  
  # un-comment when finished debugging - this will work as long as eff_layer only has one layer
  colnames(eff_df) <- c('x', 'y', 'layer')
  
  #get species list from length of community list
  species_list <- vector()
  for (i in 1:length(community)){species_list[i] <- paste0("Sp",i)}
  
  #for each species on the list, extract the relevant model outputs (this allows for some models to fail for some species)
  
  community_preds <- list()
  
  for (j in 1:length(species_list)){
    
    species <- species_list[j]
    #get model outputs for this species
    models_to_read <- grep(paste0(species,"_"), models)
    
    #load outputs into list and extract prediction table
    model_outputs <- list()
    idx <- 1
    for (k in models_to_read){
      try(load(paste0(sdm_path, models[k])))
      model_type <- model_output$model
      if (model_type != "rf"){
        model_preds <- model_output$predictions}
      if (model_type == "rf"){
        model_preds <- model_output$predictions[,names(model_output$predictions) %in% c("x", "y", "mean.1", "sd", "DECIDE_score.1")]
        names(model_preds) <- c("x", "y", "mean", "sd", "DECIDE_score")
      }
      model_preds$mean <- model_preds$mean*(1-prevalence_vec[j]) #weight prevalence by rarity
      model_preds$DECIDE_score <- model_preds$mean*model_preds$sd #recalculate with prevalence weighted by rarity
      model_outputs[[idx]] <- model_preds
      names(model_outputs)[idx] <- model_type
      idx <- idx + 1
    }
    
    #average model outputs (note - not weighted by AUC currently)
    mod_average <- Reduce(`+`, model_outputs) / length(model_outputs)
    
    #multiply mean by 1-prevalence - upweights rare species - now done per model above
    #try(mod_average$mean <- mod_average$mean*(1-prevalence_vec[j]))
    
    #store only the model average for now - could edit to store the individual model outputs if needed
    if(is.null(nrow(mod_average))){community_preds[[j]] <- NULL} else { community_preds[[j]] <- mod_average; names(community_preds)[j] <- species}
    
    
  }
  
  # print a little output message
  print(paste("Finished processing all species in community", strsplit(basename(as.character(community_file)),"\\.")[[1]][1]))
  
  #average across all species in community (which can be modelled) to obtain a single prevalence, uncertainty and DECIDE score
  
  community_preds <- Filter(length, community_preds)
  
  community_scores <- Reduce(`+`, community_preds)/length(community_preds)
  
  if (method == "none"){
    cell_weights <- eff_df$layer/sum(eff_df$layer, na.rm=TRUE)
    #assign NA values the average weight
    cell_weights[is.na(cell_weights)] <- mean(cell_weights, na.rm= TRUE)
    #sample new locations according to cell weights
    new_locs <- sample(1:nrow(eff_df), size = n, replace = FALSE, prob = cell_weights)
    new_coords <- eff_df[new_locs, 1:2]
  }
  
  if (method == "uncertainty"){    #merge with existing sampling bias if uptake isn't NULL
    if(is.null(uptake)){
      cell_weights <- community_scores$sd/sum(community_scores$sd, na.rm=TRUE)
      #assign NA values the average weight
      cell_weights[is.na(cell_weights)] <- mean(cell_weights, na.rm= TRUE)
      #sample new locations according to cell weights
      new_locs <- sample(1:nrow(community_scores), size = n, replace = FALSE, prob = cell_weights^probability_weight_adj)
      new_coords <- community_scores[new_locs, 1:2]
    }
    if(!is.null(uptake)){
      #combine effort and score dataframes
      comb_df <- merge(eff_df, community_scores, by = c("x", "y"))
      #standardise both effort and score to 0 to 1
      comb_df[,3:6] <- apply(comb_df[,3:6],2,FUN = function(x) {x/max(x, na.rm=TRUE)})
      comb_df$comb_weight <- (comb_df$layer*(1-uptake))+(comb_df$sd*uptake)
      cell_weights <- comb_df$comb_weight/sum(comb_df$comb_weight, na.rm=TRUE)
      #assign NA values the average weight
      cell_weights[is.na(cell_weights)] <- mean(cell_weights, na.rm= TRUE)
      #sample new locations according to cell weights
      new_locs <- sample(1:nrow(comb_df), size = n, replace = FALSE, prob = cell_weights)
      new_coords <- comb_df[new_locs, 1:2]
    }
  }
  
  if (method == "prevalence"){
    #merge with existing sampling bias if uptake isn't NULL
    if(is.null(uptake)){
      cell_weights <- community_scores$mean/sum(community_scores$mean, na.rm=TRUE)
      #assign NA values the average weight
      cell_weights[is.na(cell_weights)] <- mean(cell_weights, na.rm= TRUE)
      #sample new locations according to cell weights
      new_locs <- sample(1:nrow(community_scores), size = n, replace = FALSE, prob = cell_weights^probability_weight_adj)
      new_coords <- community_scores[new_locs, 1:2]
    }
    if(!is.null(uptake)){
      #combine effort and score dataframes
      comb_df <- merge(eff_df, community_scores, by = c("x", "y"))
      #standardise both effort and score to 0 to 1
      comb_df[,3:6] <- apply(comb_df[,3:6],2,FUN = function(x) {x/max(x, na.rm=TRUE)})
      comb_df$comb_weight <- (comb_df$layer*(1-uptake))+(comb_df$mean*uptake)
      cell_weights <- comb_df$comb_weight/sum(comb_df$comb_weight, na.rm=TRUE)
      #assign NA values the average weight
      cell_weights[is.na(cell_weights)] <- mean(cell_weights, na.rm= TRUE)
      #sample new locations according to cell weights
      new_locs <- sample(1:nrow(comb_df), size = n, replace = FALSE, prob = cell_weights)
      new_coords <- comb_df[new_locs, 1:2]
    }
  }
  
  if (method == "unc_plus_prev"){
    #merge with existing sampling bias if uptake isn't NULL
    if(is.null(uptake)){
      cell_weights <- community_scores$DECIDE_score/sum(community_scores$DECIDE_score, na.rm=TRUE)
      #assign NA values the average weight
      cell_weights[is.na(cell_weights)] <- mean(cell_weights, na.rm= TRUE)
      #sample new locations according to cell weights
      new_locs <- sample(1:nrow(community_scores), size = n, replace = FALSE, prob = cell_weights^probability_weight_adj)
      new_coords <- community_scores[new_locs, 1:2]
    }
    if(!is.null(uptake)){
      #combine effort and score dataframes
      comb_df <- merge(eff_df, community_scores, by = c("x", "y"))
      #standardise both effort and score to 0 to 1
      comb_df[,3:6] <- apply(comb_df[,3:6],2,FUN = function(x) {x/max(x, na.rm=TRUE)})
      comb_df$comb_weight <- (comb_df$layer*(1-uptake))+(comb_df$DECIDE_score*uptake)
      cell_weights <- comb_df$comb_weight/sum(comb_df$comb_weight, na.rm=TRUE)
      #assign NA values the average weight
      cell_weights[is.na(cell_weights)] <- mean(cell_weights, na.rm= TRUE)
      #sample new locations according to cell weights
      new_locs <- sample(1:nrow(comb_df), size = n, replace = FALSE, prob = cell_weights)
      new_coords <- comb_df[new_locs, 1:2]
    }
  }
  
  if (method == "coverage"){
    # if(is.null(eff_layer)) stop("Cannot adaptively sample using method 'coverage' without specifying an existing sampling raster (i.e. currently eff_layer = NULL)")
    
    require(tidyverse)
    ## create the observations dataframe to find unsampled cells
    # get all the records that have been made in a community
    # this is equivalent to eff_df but for the simulated data rather than the real butterfly data
    obvs <- na.omit(do.call(rbind, lapply(community, function(x) x$observations)))
    
    # get the number of records in each recorded grid cell
    obvs_recs <- obvs %>% mutate(x=lon,y=lat,lon=NULL,lat=NULL) %>% 
      group_by(x, y) %>% tally
    
    #combine effort and score dataframes to get the empty cells as well
    comb_df <- merge(community_scores, obvs_recs, by = c("x", "y"), all.x = TRUE, sort = FALSE)
    comb_df$n[is.na(comb_df$n)] <- 0 # define NAs as 0; NAs are where there are no records
    
    empty_cells <- comb_df[comb_df$n == 0,]
    
    if(is.null(uptake)){
      new_locs <- sample(1:nrow(empty_cells), size = n, replace = FALSE)#sample from empty cells with equal prob
      new_coords <- empty_cells[new_locs,1:2]
    }
    if(!is.null(uptake)){
      new_locs_empty <- sample(1:nrow(empty_cells), size = n*uptake, replace = FALSE)#sample from empty cells with equal prob
      new_coords_empty <- empty_cells[new_locs_empty,1:2]
      ## sample from background effort - in pattern dictated by underlying sampling
      # calculate background sampling probability - copied from 'none'
      cell_weights <- eff_df$layer/sum(eff_df$layer, na.rm=TRUE)
      #assign NA values the average weight
      cell_weights[is.na(cell_weights)] <- mean(cell_weights, na.rm= TRUE)
      #sample new locations according to cell weights
      new_locs_effort <- sample(1:nrow(eff_df), size = n*(1-uptake), replace = FALSE, prob = cell_weights)
      new_coords_effort <- eff_df[new_locs_effort, 1:2]
      new_coords <- rbind(new_coords_empty, new_coords_effort) # bind the new coords together
    }
  }
  
  if (method == "unc_plus_recs"){
    
    require(tidyverse)
    
    # get all the records that have been made in a community
    # this is equivalent to eff_df but for the simulated data rather than the real butterfly data
    obvs <- na.omit(do.call(rbind, lapply(community, function(x) x$observations)))
    
    # get the number of records in each recorded grid cell
    num_recs <- obvs %>% mutate(x=lon,y=lat,lon=NULL,lat=NULL) %>% 
      group_by(x, y) %>% tally 
    
    #combine effort and score dataframes
    comb_df <- merge(community_scores, num_recs, by = c("x", "y"), all.x = TRUE, sort = FALSE)
    comb_df$n[is.na(comb_df$n)] <- 0 # define NAs as 0; NAs are where there are no records
    comb_df$n <- comb_df$n+1 # +1 to all so that when dividing by areas with no records, just get the unaltered uncertainty score
    
    # weighting (current score in DECIDE) is done by unc * 1 over number of records in a cell
    # well actually, it's 1/time since last record - but that's not possible in our framework unless we use the raw 
    # data which might cause problems...
    comb_df$unc_recs <- comb_df$sd*(1/comb_df$n)
    
    if(is.null(uptake)){
      cell_weights <- comb_df$unc_recs/sum(comb_df$unc_recs, na.rm=TRUE)
      #assign NA values the average weight
      cell_weights[is.na(cell_weights)] <- mean(cell_weights, na.rm= TRUE)
      #sample new locations according to cell weights
      new_locs <- sample(1:nrow(comb_df), size = n, replace = FALSE, prob = cell_weights^probability_weight_adj)
      new_coords <- comb_df[new_locs, 1:2]
    } else if(!is.null(uptake)) {
      # combine with the specified effort layer dataset above so can see what happens
      # if people carry on with business as usual instead of doing the sampling
      comb_df_eff <- merge(eff_df, comb_df, by = c("x", "y"), all.x = TRUE)
      comb_df_eff[,3:8] <- apply(comb_df_eff[,3:8],2,FUN = function(x) {x/max(x, na.rm=TRUE)})
      comb_df_eff$comb_weight <- (comb_df_eff$layer*(1-uptake))+(comb_df_eff$unc_recs*uptake) # change the influence of adaptive sampling
      cell_weights <- comb_df_eff$comb_weight/sum(comb_df_eff$comb_weight, na.rm=TRUE)
      #assign NA values the average weight
      cell_weights[is.na(cell_weights)] <- mean(cell_weights, na.rm= TRUE)
      #sample new locations according to cell weights
      new_locs <- sample(1:nrow(comb_df_eff), size = n, replace = FALSE, prob = cell_weights)
      new_coords <- comb_df_eff[new_locs, 1:2]
    }
    
  }
  
  
  
  ##sample new observations from coordinates
  
  
  new.obs <- list()
  for (i in 1:length(community)){
    # print(paste("NAs in coordinates? =", as.character(any(is.na(new_coords)))))
    observations <- data.frame(sp::coordinates(new_coords)[,1:2])
    observations$Real <- raster::extract(community[[i]]$pres_abs, observations)
    observations <- observations[observations$Real == 1,]#remove absences to create presence-only data?
    observations$Observed <- observations$Real * (rbinom(nrow(observations),1,0.2)) # are species detected?
    observations <- observations[!is.na(observations$x),]#remove missing locs
    if(nrow(observations)>0) observations[observations == 0] <- NA # set places where no individual was observed to NA - $Observed added to dea
    names(observations) <- c("lon", "lat", "Real", "Observed")
    new.obs[[i]] <- list()
    new.obs[[i]]$observations <- observations
  }
  
  community_AS <- lapply(seq_along(community), function(x) list(observations = rbind(community[[x]]$observations, new.obs[[x]]$observations), variables = community[[x]]$variables, model_variables = community[[x]]$model_variables, AS_method = method))
  
  community_name <- strsplit(basename(as.character(community_file)),"\\.")[[1]][1]
  
  # # This is to create a separate location for adaptively sampled data if we want - probably not worth it.
  # dir.create(paste0(outPath, version_name, "adaptive_sampling_data/"), recursive = TRUE)
  
  saveRDS(community_AS, file = paste0(outPath, AS_version, '_', community_name, "_AS_", method, ".rds"))
  
}