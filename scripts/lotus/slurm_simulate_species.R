simulate_species <- function(env_data, extent = NULL, n = 10, outPath, seed = NULL, n_env = NULL, beta = 0.5, alpha = -0.05, max_samp = 1000, det_prob = 0.5, effort = NULL, weight_adj = 1, background = NULL){
  
  library(raster)
  library(virtualspecies)
  
  #read in data
  env <- raster::stack(env_data)
  
  #set seed if specified
  if(!is.null(seed)){set.seed(seed)}
  
  #crop to extent if specified
  if(!is.null(extent)){
    e <- as(extent, "SpatialPolygons")
    sp::proj4string(e) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    
    e.geo <- sp::spTransform(e, CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs"))
    
    env_extent <- raster::crop(env, e.geo)
    
  } else {env_extent <- env}
  
  #subset env data layers
  
  #decide on number of layers to use or range of layers from which to select, if no limit given use all layers in raster
  if(is.null(n_env)){env_min = nlayers(env_extent); env_max = nlayers(env_extent)} else if (length(n_env) == 1) {env_min = n_env; env_max = n_env} else if(length(n_env) == 2) {env_min = n_env[1]; env_max = n_env[2]} else if (length(n_env) > 2) {stop("n_env greater than 2: Input either a single number or a range from which to select the number of environmental layers to use in species generation")}
  
  #set background if given, can indicate a layer in env_data or be a filepath to a raster
  if(is.numeric(background)){
    bg_layer <- env_extent[[background]]
  } else if (is.character(background) & !grepl("\\.", background)) {bg_layer <- raster::subset(env_extent, background)} else if (is.character(background) & grepl("\\.", background)) {bg_layer <- raster::raster(background)} else {bg_layer <- NULL}
  
  #extract effort layer from raster if provided (note currently uses layers in existing raster stack, could read in other layers)
  if(is.numeric(effort)){eff_layer <- env_extent[[effort]]} else if(is.character(effort) & !grepl("\\.", effort)) {eff_layer <- raster::subset(env_extent,effort)} else if (is.character(effort) & grepl("\\.", effort)) {eff_layer <- raster::raster(effort)} else  {eff_layer <- NULL}
  
  if(is.null(eff_layer)){eff_weights <- (env_extent[[1]]*0)+1} else if (is.null(bg_layer)){
    eff_weights <- eff_layer/weight_adj} else {eff_weights <- (bg_layer/bg_layer) + (eff_layer/weight_adj)}
  
  
  community <- list()
  
  #for each species generate observations
  for (i in 1:n){
    #subset env raster
    my.stack <- env_extent[[sample(1:nlayers(env_extent),size = runif(1,env_min,env_max), replace = FALSE)]]
    #generate a suitability raster
    my.pca.species <- generateSpFromPCA(raster.stack = my.stack, sample.points=TRUE, nb.points = 10000, plot = FALSE)
    #convert to presence-absence
    pa <- convertToPA(my.pca.species, beta = beta, alpha = alpha, plot = FALSE)
    #extract prevalence
    prevalence <- as.numeric(pa$PA.conversion[5])
    #determine maximum number of observations based on prevalence
    #max_obs <- round(prevalence*max_samp)
    max_obs <- max_samp #set max no of observations - could use data?
    #sample observations based on bias and detection prob
    occs <- sampleOccurrences(pa, n = max_obs, type = "presence-absence", detection.probability = det_prob, bias = "manual", weights = eff_weights)
    #rename columns of occurrences data
    if(nrow(occs$sample.points) > 0){names(occs$sample.points) <- c("lon", "lat", "Real", "Observed")}
    #subset to PO data
    occs$sample.points <- occs$sample.points[occs$sample.points$Real == 1,]
    occs$sample.points$Observed[occs$sample.points$Observed == 0] <- NA
    #store required outputs to list
    community[[i]] <- list(true_prob_occ = pa$probability.of.occurrence, pres_abs = pa$pa.raster, observations = occs$sample.points, variables = pa$details$variables, prevalence = prevalence)
  }
  
  #return(community)
  
  community_name <- paste0("community_",seed,"_", n, "_sim")
  
  if(!dir.exists(paste0(outPath, community_name,"/"))){
    dir.create(paste0(outPath, community_name,"/"))
  }
  
  saveRDS(community, file = paste0(outPath,community_name,"/", community_name, ".rds"))
}

library(rslurm)

dirs <- config::get("LOTUSpaths_sim")

n_communities = 50

pars <- data.frame(env_data = "envdata_1km_no_corr_noNA.grd",outPath = dirs$outpath, seed = 1:n_communities, max_samp = 10000, n_env = 10, n = 50, effort = "butterfly_1km_effort_layer.grd", background = "MeanDiRange")

sjob <- slurm_apply(simulate_species, pars, 
                    jobname = 'sim_spp',
                    nodes = nrow(pars), 
                    cpus_per_node = 1, 
                    submit = TRUE,
                    slurm_options = list(partition = "test",
                                         time = "0:59:59",
                                         mem = "10000",
                                         output = "sim_spp_%a.out",
                                         error = "sim_spp_%a.err"),
                    sh_template = "jasmin_submit_sh.txt")

