
####    Processing environmental data    ####

library(raster)
library(sf)


# location to store joint files
out_dir <-'P:/07543_DECIDE/Data/WP1 modelling/'


### cropping raster to GB

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

####   getting 1km raster    ####

# load in lcm target class, 1km
lcm_1k <- raster::stack('data/environmental_data/1km_raster/data/LCM2015_GB_1km_percent_cover_target_class.tif')
names(lcm_1k) <- c('Broad_wood', 'conif_wood', 'arable', 'impr_grass', 'neutr_grass', 'calc_grass', 'acid_grass',
                   'marsh', 'heather', 'heather_grass', 'bog', 'inland_rock', 'saltwater', 'freshwater',
                   'sup_lit_rock', 'sup_lit_sed', 'lit_rock', 'lit_sed', 'saltmarsh', 'urban', 'suburban')

names(lcm_1k) <- names(envdat)[1:21]

# mask to only GB
lcm_1k_cropped <- mask(lcm_1k, gb_map)

# get elevation at 1km
elev_100 <- raster("C:/Users/thoval/OneDrive - UKCEH/Documents/DECIDE/DECIDE_WP1/Data/raw_data/environmental/Copernicus_Elevation/elevation_UK.tif")
elev_1km <- raster::aggregate(elev_100, fact = 10, fun = 'mean')
elev_1km

beginCluster()
elev <- raster::projectRaster(elev_1km, lcm_1k, method = 'bilinear')
endCluster()

elev_1km_GB <- mask(elev, gb_map)
plot(elev_1km_GB)


## weather variables
haduk_100 <- envdat_cropped[[22:30]]
haduk_1k <- raster::aggregate(haduk_100, fact=10, FUN='mean')
