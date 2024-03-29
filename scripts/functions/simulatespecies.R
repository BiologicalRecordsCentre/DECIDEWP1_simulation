simulate_species <- function(env_data, extent = NULL, n = 10, outPath, seed = NULL, n_env = NULL, beta = 0.5, alpha = -0.05, max_samp = 1000, det_prob = 0.5, effort = NULL, weight_adj = 10){
  
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
  
  #extract effort layer from raster if provided (note currently uses layers in existing raster stack, could read in other layers)
  if(is.numeric(effort)){eff_layer <- env_extent[[effort]]} else if(is.character(effort)) {eff_layer <- subset(env_extent,effort)} else {eff_layer <- NULL}
  
  if(is.null(eff_layer)){eff_weights <- (env_extent[[1]]*0)+1} else {
    eff_weights <- eff_layer/weight_adj}
  
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
    max_obs <- round(prevalence*max_samp)
    #sample observations based on bias and detection prob
    occs <- sampleOccurrences(pa, n = max_obs, type = "presence only", detection.probability = det_prob, bias = "manual", weights = eff_weights)
    #rename columns of occurrences data
    if(nrow(occs$sample.points) > 0){names(occs$sample.points) <- c("lon", "lat", "Real", "Observed")}
    #store required outputs to list
    community[[i]] <- list(true_prob_occ = pa$probability.of.occurrence, pres_abs = pa$pa.raster, observations = occs$sample.points, variables = pa$details$variables)
  }
  
  return(community)
}