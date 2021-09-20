# calculating number of visits to unique locations in the butterfly dataset.

library(tidyverse)


# function to floor the coordinates rather than round
signif.floor <- function(x, n){
  pow <- floor( log10( abs(x) ) ) + 1 - n
  y <- floor(x / 10 ^ pow) * 10^pow
  # handle the x = 0 case
  y[x==0] <- 0
  y
}

b <- read_csv('../DECIDE_WP1/data/edited_insect_data/butterfly/butterfly_EastNorths_no_duplicates.csv')
head(b)
head(b$lon)

# number of records per species per year
b %>% 
  mutate(rlon = signif.floor(lon, 3), # three is 1km
         rlat = signif.floor(lat, 3)) %>% 
  distinct(sp_n, rlon, rlat, yr) %>% 
  group_by(sp_n, yr) %>% 
  tally %>% 
  ggplot(aes(x = reorder(sp_n, n), y = n)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1))


## number of unique locations
b %>% 
  mutate(rlon = signif.floor(lon, 3), # three is 1km
         rlat = signif.floor(lat, 3)) %>% 
  distinct(rlon, rlat) %>% # number of distinct sites at 1km
  dim()
# 17297


## number of unique locations visited in a single year
nu_yr <- b %>% 
  mutate(rlon = signif.floor(lon, 3), # three is 1km
         rlat = signif.floor(lat, 3)) %>% 
  distinct(rlon, rlat, yr) %>% # distinct sites in each year
  mutate(loc_id = paste(rlon, rlat, sep = '_')) %>% 
  group_by(yr) %>% 
  summarise(visits = length(unique(loc_id)))
head(nu_yr)

mean(nu_yr$visits[nu_yr$yr!=2016])


## test the raster
library(raster)

t <- raster::stack(list.files("C:/Users/thoval/OneDrive - UKCEH/Documents/DECIDE/DECIDE_WP1/data/raw_data/environmental/HadUK_dat/", full.names = TRUE)[1])[[1]]
plot(t)


## crop to GB
# download map GB
uk_map <- st_as_sf(getData("GADM", country = "GBR", level = 1, path='data/environmental_data'))
uk_map <- st_transform(uk_map, 27700)

# remove nrothern ireland
gb_map <- uk_map[uk_map$NAME_1 != 'Northern Ireland',]

# check
plot(st_geometry(gb_map))

# convert to spatial for use in raster::mask()
gb_mask <- as_Spatial(gb_map)
gb_mask

# mask elevation
m_gb <- raster::mask(t, gb_mask[1])
plot(m_gb)

as.data.frame(m_gb) %>% na.omit %>% 
  dim()

m_gb
