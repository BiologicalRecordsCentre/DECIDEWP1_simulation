

library(raster)
library(tidyverse)

source('scripts/functions/simulatespecies.R')

### step 1 load data

# read in data
epath <- 'P:/07543_DECIDE/Data/WP1 modelling/'
env_data <- raster::stack(paste0(epath, "edat_nocorrs_nosea.gri"))[[19:23]]
names(env_data)
env_dat <- env_data[[c(19:22,31)]]


# crop to a region
ext_h <- extent(matrix(c(-3,53.5, -2,54.5), ncol = 2))
e <- as(ext_h, "SpatialPolygons")
sp::proj4string(e) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
e.geo <- sp::spTransform(e, CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs"))
env_data_cropped <- raster::crop(env_dat, e.geo)
plot(env_data_cropped)


### step 2 simulate species
sim_spp <- simulate_species(env_data = env_data_cropped,
                            n=2,
                            effort = 'suburban')

### step 3 convert to cpa() format

t <- sim_spp[[2]]
