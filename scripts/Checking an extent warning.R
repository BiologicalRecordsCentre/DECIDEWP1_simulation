
## checking warning error that was produced whenever using the butterfly effort layer
## the warning was rasters of different extents - returning overlap. 
## judging from below, this isn't a problem, but wanted to check.

library(raster)

b_layer <- raster::stack('data/environmental_data/butterfly_1km_effort_layer.grd')
plot(b_layer)

e_layer <- raster::stack('data/environmental_data/envdata_1km_no_corr_noNA.grd')
plot(e_layer[[1]])

cropped_e <- raster::crop(e_layer, b_layer)

par(mfrow=c(1,2))
plot(e_layer[[1]])
plot(cropped_e[[1]])