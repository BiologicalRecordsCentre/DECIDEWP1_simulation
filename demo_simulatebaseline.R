#' # Example of basic virtual species simulation on subset of GB focused on northern England.
#' 
#' The steps are as follows:
#' 
#' 1. Define area of interest and extract environmental rasters from set used in DECIDE modelling for this area. 
#' 
#' 2. Subset environmental raster layers to smaller set for this example
#' 
#' 3. Use virtualspecies package to generate a virtual species distribution
#' 
#' 4. From the species distribution determine where the species is present and where it is absent (i.e. create binary map)
#' 
#' 5. Sample the virtual species distribution according to existing sampling effort pattern (TO EXTRACT - currently uses suburban area to mimic variation in effort)



#' ## 1. Define area of interest and extract environmental rasters

library(raster)
library(virtualspecies)

set.seed(1000)#set seed for test runs to check performance

#set path to DECIDE env data (mapped mine to F as issues reading off Wallingford P drive for some reason)
epath <- "F:Data\\WP1 modelling\\"
epath <- 'P:/07543_DECIDE/Data/WP1 modelling/'

# read in data
env_data <- raster::stack(paste0(epath, "edat_nocorrs_nosea.gri"))

# crop to region used by Thomas in model testing using Thomas' code
ext_h <- extent(matrix(c(-4,53, 0.2,54.5), ncol = 2))
e <- as(ext_h, "SpatialPolygons")
sp::proj4string(e) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

e.geo <- sp::spTransform(e, CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs"))

hbv_y <- raster::crop(env_data, e.geo)

par(mfrow=c(1,1))
plot(hbv_y[[1]], main = names(hbv_y[[1]]))


#' ## 2. Subset environmental layers

#' Two options given here - either randomly select a random number of layers or use all layers (only useful for testing?)

#random selection of covariates to test in virtual species generation - randomly samples between 5 and 25 layers of the raster stack (i.e. not all species would be based on same set of covariates, or even the same number of layers)
my.stack <- hbv_y[[sample(1:nlayers(hbv_y),size = runif(1,5,25), replace = FALSE)]]

#use all PCA layers - test to see how long this takes with PCA generation. Not too slow with all layers on regional subset = no risk to allowing quite a few layers in random generation
#my.stack = hbv_y

#' We can extract a summary of the layers randomly selected 
my.stack

#' ## 3. Use virtualspecies package to generate virtual species

#' We can either generate species with random properties (i.e. the package will pick for us where the species lives in environmental space) or provide more specific instructions

#' Firstly generate a species randomly with no custom input

#generate a virtual species using the PCA approach - because of high resolution use the sample.points argument to pick random sample of grid cells to use in the PCA (subsetted areas has 4.7 million cells, pick 10,000 to do PCA for speed)
my.pca.species <- generateSpFromPCA(raster.stack = my.stack, sample.points=TRUE, nb.points = 10000)

#extract summary of virtual species
my.pca.species

#' If required, we could pick species with certain niche breadths, the package allows either broad or narrow niche widths to be specified (otherwise the package will select randomly, in the example above we randomly generated a generalist species

narrow.species <- generateSpFromPCA(raster.stack = my.stack, sample.points=TRUE, nb.points = 10000, niche.breadth = "narrow")

narrow.species

#' This species has a very narrow niche breadth and is restricted to a small subset of the environmental space
#' 
#' It is possible to customise very specifically the environmental characteristics, but when simulating many species it may be more reasonable to let the package choose for us
#' 
#' To demonstrate, lets simulate 10 species using randomly selected niche breadths and see if the range of distributions seems reasonable
#' 

virt_comm1 <- list()
for (i in 1:10){
my.stack <- hbv_y[[sample(1:nlayers(hbv_y),size = runif(1,5,25), replace = FALSE)]]
my.pca.species <- generateSpFromPCA(raster.stack = my.stack, sample.points=TRUE, nb.points = 10000)
virt_comm1[[i]] <- my.pca.species
}

par(mfrow=c(5,2))
par(mar= c(1,1,1,1))
for(i in 1:10){
  plot(virt_comm1[[i]])
}

#' Seems like a reasonable mix of species with wide and narrow niche breadths
#'
#' Alternatively we can generate species without the PCA step by generating random responses to each of the environmental gradients (functions are available to parameterise these responses using e.g. known species niches but this is probably a level of complexity above what we need in this project).
#' 
virt_comm2 <- list()
for (i in 1:10){
my.stack <- hbv_y[[sample(1:nlayers(hbv_y),size = runif(1,5,25), replace = FALSE)]]
random.sp <- generateRandomSp(my.stack, approach = "response", convert.to.PA = FALSE, realistic.sp = TRUE, plot = FALSE)
virt_comm2[[i]] <- random.sp
}

par(mfrow=c(5,2))
par(mar= c(1,1,1,1))
for(i in 1:10){
  plot(virt_comm2[[i]])
}

#' Seems to struggle to come up with reasonable species distributions, PCA approach seems sensible

#' 
#' 
#' 4. Create binary map
#' 
#' Once the underlying distribution/environmental suitability is defined we can then determine where the species is present (i.e. the species won't be present in all areas where suitability is high and may be present in some areas where suitability is only moderate)
#' 
#'  We determine where the species is present probabilistically so that the regions of highest suitability are those with highest probability of a species being present. By default a logistic conversion is used to convert environmental suitability to probability of occurrence. The shape of the logistic relationship means that the differences between suitability scores have most impact on probability of occurrence at medium values. This conversion is customisable but the default will randomly assign conversion parameters.
#'  
#'  We can see what this looks like for our virtual community using randomly assigned conversions
#' 
par(mfrow=c(5,2)) 
pa <- list()
for(i in 1:10){
pa[[i]] <- convertToPA(virt_comm1[[i]], plot = FALSE)
plot(pa[[i]]$pa.raster)
}

#' Alternatively it may be more sensible to use the same conversion for each species to ensure environmental suitability is converted in a standardised way 
#' 
par(mfrow=c(5,2)) 
pa <- list()
for(i in 1:10){
  pa[[i]] <- convertToPA(virt_comm1[[i]], plot = FALSE, beta = 0.5, alpha = -0.05)
  plot(pa[[i]]$pa.raster)
}

#' It is also possible to modify this conversion to reach a specific level of coverage/prevalence i.e. to simulate a species that covers 5% or 50% of the available area. However, by using the standard conversion we will already achieve species with different prevalences so it may not be neccessary to adjust this. 
#' 
#' 
#' 5. Sample occurrences
#' 
#' Let's look at the first of our species again
#' 
plot(pa[[1]])
#'    
#' We now want to simulate sampling of these presences in a way that mirrors real life sampling:
#' 
#' 1. Sampling is biased to certain regions (easily accessible, high human population density etc)
#' 
#' 2. Detection probability is less than 1 (note that initially we'll assume detection probability is the same a) between species b) across space - these are big assumptions that we'll come back to later and assess)
#' 
#' 3. We sample presences only, but use presences from other species to create pseudo-absences
#' 
#' The virtualspecies package requires a maximum number of observations to sample which are then thinned by removing absences. Note that an alternative method to generate these presences would be to use the probability of occurrence or suitability rasters to create an intensity raster (effectively interpreted as the density of individuals per unit area). We could then use point process theory to generate presences as realisations of an (inhomogenous or log Gaussian Cox) point process. Point process theory states that the higher the intensity, the more points (individuals) per unit area are expected. However, sampling using this method is likely to be computationally inefficient with a high resolution problem so for now we'll ignore this (but note that there are advantages of a point process framework we might want to come back to later, not least the potential to combine presence-only, presence-absence data and abundance in continuous space and the ability to generate a random number of presences conditional on the intensity).
#' 
#' Ok, so let's simulate some presence-only data derived from species 1, assuming a detection probability of 0.5.
#' 
#' 
par(mfrow=c(1,1))
sp.obs1 <- sampleOccurrences(pa[[1]], n = 100, type = "presence only", detection.probability = 0.5)

length(sp.obs1$sample.points$Observed[!is.na(sp.obs1$sample.points$Observed)])

#' Most observations in the east where suitability is higher, note that due to the detection probability being less than 1 the number of observations retrieved is less than 100
#' 
#' However, in the current set up the number of observations made is deterministic and not a function of prevalence (for example we might expect a widespread species to have more observations than a species restricted to moorlands). In addition, there is currently no sampling bias, and equal effort is expended across the region. 
#' 
#' Next step - unequal effort
#' 
#' Firstly, we can introduce a sampling bias related to environment - for example we could assume that high suburban sites are more commonly visited than areas with low surburan area
#' 
suburban <- hbv_y[[21]]

#we can set weights so that the highest suburban areas are 10 times more likely to be sampled than the lowest suburban areas (on the original scale suburban areas were 100 times more likely to be sampled which may be extreme)

sub_weight <- suburban/10

#' Now we can use this weight to alter our sampling 
#' 
sp.obs2 <- sampleOccurrences(pa[[1]], n = 100, type = "presence only", detection.probability = 0.5, bias = "manual", weights = sub_weight)
#'  
#' Real life sampling biases are often really challenging to explain as they are complex combinations of the many motivations of observers. However, we can approximate these by using known sampling patterns from existing recorders to simulate realistic patterns in sampling effort
#' 
#' One other element to add would be to specify n (the number of observations prior to thinning due to detection probability) as proportional to prevalence i.e. those species that are more widespread would have a higher potential number of observations. This is currently not included as part of the virtualspecies package but would be easy to implement    
#' 
#extract prevalence for all species in virt_comm1
prev_vec <- vector()
for (i in 1:length(pa)){
  prev_vec[i] <- as.numeric(pa[[i]]$PA.conversion[5])
}
prev_vec
#' Prevalence ranges between 0.021 and 0.977 (possibly a bit high overall? reasonable for trials)
#' 
#' Sample proportional to prevalence, weighted by suburban area
#' 
par(mfrow=c(5,2))
sp.obs <- list()
for (i in 1:10){
max_obs <- round(prev_vec[i]*1000)
sp.obs[[i]] <- sampleOccurrences(pa[[i]], n = max_obs, type = "presence only", detection.probability = 0.5, bias = "manual", weights = sub_weight)
names(sp.obs[[i]]$sample.points) <- c("lon", "lat", "Real", "Observed")
}

save.image(file = "virt_comm_10spp.Rdata")
