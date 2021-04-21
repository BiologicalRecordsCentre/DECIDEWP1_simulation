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
#' 4. Sample the virtual species distribution according to existing sampling effort pattern (TO EXTRACT)



#' ## 1. Define area of interest and extract environmental rasters

library(raster)
library(virtualspecies)

set.seed(1000)#set seed for test runs to check performance

#set path to DECIDE env data (mapped mine to F as issues reading off Wallingford P drive for some reason)
epath <- "F:Data\\WP1 modelling\\"

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
#' 
#' 4. Sample from species distribution
#' 
#' Once the underlying distribution/environmental suitability is defined we can then sample from it to determine where the species is present and observed
#' 
#'  We do this probabilistically so that the regions of highest suitability are those with highest likelihood of a species being present. By default a logistic conversion is used to convert environmental suitability to probability of occurrence. The shape of the logistic relationship means that the differences between suitability scores have most impact on probability of occurrence at medium values. This conversion is customisable but the default will randomly assign conversion parameters.
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

#'    
#'    
#'    