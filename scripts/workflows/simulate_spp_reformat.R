

library(raster)
library(tidyverse)
library(dismo)
library(randomForest)

source('scripts/functions/simulatespecies.R')
source('scripts/functions/reformat_simulated_data.R')
source('scripts/functions/Edited_Rob_Functions.R')

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
system.time(sim_spp <- simulate_species(env_data = env_data_cropped,
                                        n=5,
                                        effort = 'suburban'))


### step 3 convert to cpa() format using new function
presences_df <- reformat_data(sim_spp, year = 2015, species_name = 'Sp')
head(presences_df)


### step 4 run presence absence function
# get unique list of names
species_list <- unique(presences_df$species)

pres_abs <- vector('list', length = length(species_list))

for(s in 1:length(species_list)){
  
  pres_abs[[s]] <- cpa(spdat = presences_df, species = species_list[s], 
                       matchPres = FALSE, nAbs = 10000,
                       minYear = 2000, maxYear = 2017, recThresh = 5,
                       screenRaster = env_data_cropped)
  
}

names(pres_abs) <- species_list


### step 5 run model
sdm_lr <- fsdm(species = species_list[1], model = "rf",
               climDat = env_data_cropped, spData = pres_abs, knots_gam = 4,
               k = 4, 
               write =  F, outPath = "C:/Users/thoval/Documents/Analyses/lr_outs/")

preds <- get_predictions(sdm_lr, 'rf', env_data_cropped)
