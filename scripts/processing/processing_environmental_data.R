
####    Processing environmental data    ####

library(raster)
library(sf)

# install.packages("ncdf4")

# location to store joint files
out_dir <-'P:/07543_DECIDE/Data/WP1 modelling/'


####   FIRST, crop 100m raster to GB    ####

# read in current file
envdat <- raster::stack(paste0(out_dir, 'edat_nocorrs_nosea.grd'))

# download map of GB
uk_map <- st_as_sf(getData("GADM", country = "GBR", level = 1, path='data/environmental_data/'))
uk_map <- st_transform(uk_map, 27700)
gb_map <- as_Spatial(uk_map[uk_map$NAME_1 != 'Northern Ireland',])
plot(gb_map)

envdat_cropped <- mask(envdat, gb_map)
envdat_cropped
plot(envdat_cropped[[31]])

# # save raster to shared drive
# writeRaster(x = envdat_cropped, 
#             filename = paste0(out_dir, "edat_nocorrs_nosea_cropped.grd"),
#             format = 'raster', overwrite = T)


envdat_cropped <- raster::stack(paste0(out_dir, "edat_nocorrs_nosea_cropped.grd"))



####   SECOND, getting a 1km raster    ####

# load in lcm target class, 1km, which was downloaded from CEH
# loaded in from shared directory
lcm_1k <- raster::stack(paste0(out_dir,'data/environmental_data/1km_raster/data/LCM2015_GB_1km_percent_cover_target_class.tif'))
# names(lcm_1k) <- c('broad_wood', 'conif_wood', 'arable', 'impr_grass', 'neutr_grass', 'calc_grass', 'acid_grass',
#                    'fen_marsh_swamp', 'heather', 'heather_grass', 'bog', 'inland_rock', 'saltwater', 'freshwater',
#                    'sup_lit_rock', 'sup_lit_sed', 'lit_rock', 'lit_sed', 'saltmarsh', 'urban', 'suburban')

names(lcm_1k) <- names(envdat_cropped)[1:21]

# mask to only GB
lcm_1k_cropped <- mask(lcm_1k, gb_map)

# get elevation at 1km
# read in from shared folder
elev_100 <- raster(paste0(out_dir,'data/environmental_data/raw_data/copernicus_elevation/elevation_UK.tif'))
elev_1km <- raster::aggregate(elev_100, fact = 10, fun = 'mean')
elev_1km

beginCluster()
elev <- raster::projectRaster(elev_1km, lcm_1k, method = 'bilinear')
endCluster()

elev_1km_GB <- mask(elev, gb_map)
plot(elev_1km_GB)


####    weather variables    ####
# to create the bioclimatic variables I have used Had UK climate data
# I have transferred these to the shared folder
list.files(paste0(out_dir, 'data/environmental_data/raw_data/HadUK_dat/'))


# now want to get thede files into a format to create the bioclim variables
path <- paste0(out_dir, 'data/environmental_data/raw_data/HadUK_dat/')

vars <- c("rainfall", "tasmax", "tasmin")

out_var <- list()

for(i in vars) {
  
  # list all the files with a given structure
  files <- list.files(pattern = i, path = path, full.names = T)
  
  # initialise list
  st_l <- list()
  
  for(j in 1:length(files)){
    
    # get the files for each year
    # need to use stack so that each month is a different layer
    env_d <- raster::stack(files[j])
    
    # rename layers 1:12 - months of the year
    names(env_d) <- paste("M", seq(1:12), sep = "_")  
    
    # store as an object to use as an index to average across later
    ind <- names(env_d)
    
    # store the object in the list to use later
    st_l[[j]] <- env_d
    
    
  }
  
  print(i)
  
  # now want to average the same layer across multiple raster stacks
  # basically because I want to get the average min, max temp and rainfall 
  # for each month over 6 years, 2010 to 2015 inclusive
  st_o <- stackApply(stack(st_l), # creates a stack of all objects in the list
                     fun = mean, # takes the mean
                     indices = ind, # provide an index to get the names
                     na.rm = T)
  out_var[[i]] <- st_o
  
}

# convert to the 19 bioclimatic variables
had_bv <- dismo::biovars(prec = out_var$rainfall, tmin = out_var$tasmin, tmax = out_var$tasmax)

# name real names
names(had_bv) <- c("AnnualTemp",
                   "MeanDiRange",
                   "Isotherm",
                   "TempSeasonality",
                   "MaxTempWarmestMonth",
                   "MinTempColdestMonth",
                   "TempAnnualRange",
                   "MeanTempWetQuarter",
                   "MeanTempDriestQuarter",
                   "MeanTempWarmQuarter",
                   "MeanTempColdQuarter",
                   "AnnualPrecip",
                   "PrecipWetMonth",
                   "PrecipDriestMonth",
                   "PrecipSeasonality",
                   "PrecipWettestQuarter",
                   "PrecipDriestQuarter",
                   "PrecipWarmQuarter",
                   "PrecipColdQuarter")

plot(had_bv[[1]])

had_bv_cropped <- mask(had_bv, gb_map)

plot(had_bv_cropped[[1]])

beginCluster()
had_bv_cropped <- raster::projectRaster(had_bv_cropped, lcm_1k_cropped, method = 'bilinear')
endCluster()


# # store the 1km UK
# writeRaster(had_bv_cropped,
#             filename = paste0(out_dir, 'data/environmental_data/had_bv_1km_national_grid.grd'),
#             format = "raster", overwrite = T)


# stack them all together
gb_1km <- stack(lcm_1k_cropped, had_bv_cropped, elev_1km_GB)
gb_1km

# get slope and aspect from the elevtion layer
slope_asp <- terrain(gb_1km[[41]], opt = c('slope', 'aspect'), unit = 'degrees', neighbors = 8)

gb_1km <- stack(gb_1km, slope_asp)

# get rid of variables that are correlated based on 100m raster 
names(gb_1km)[names(gb_1km) %in% names(envdat)]

gb_1km_nocorr <- gb_1km[[c(names(gb_1km)[names(gb_1km) %in% names(envdat)])]]
gb_1km_nocorr

# # save the new 1km raster
# writeRaster(gb_1km_nocorr,
#             filename = paste0(out_dir, 'data/environmental_data/envdata_1km_no_corr.grd'),
#             format = "raster", overwrite = T)
