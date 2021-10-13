
## produce rank abundance curves for moths and butterflies
# data not held in this directory

library(tidyverse)

# butterflies
bdf <- read_csv("~/DECIDE/DECIDE_WP1/data/edited_insect_data/butterfly/butterfly_EastNorths_no_duplicates.csv") %>% 
  as.data.frame()
head(bdf)
bdf$lon[1]
round(bdf$lon[1], -2)

# get total number of visited locations
tot_locs <- (bdf %>% 
  mutate(rlon = round(lon, -2), 
         rlat = round(lat, -2)) %>% 
  dplyr::select(rlon, rlat) %>% 
  distinct() %>% 
  dim())[1]

# number of unique locations per species
bs <- bdf %>% 
  mutate(rlon = round(lon, -2), 
         rlat = round(lat, -2)) %>% 
  dplyr::select(rlon, rlat, sp_n) %>% 
  distinct(.keep_all = F)

# calculate prevalence for each species
b_sum <- bs %>% 
  group_by(sp_n) %>% 
  summarise(un_locs = length(rlon),
            prev = un_locs/tot_locs) %>% 
  ungroup()
b_sum

ggplot(b_sum, aes(x = reorder(sp_n, -prev), y = prev)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle('Butterfly')


# moths
mdf <- read_csv("~/DECIDE/DECIDE_WP1/data/edited_insect_data/moth/DayFlyingMoths_EastNorths_no_duplicates.csv") %>% 
  as.data.frame()
head(mdf)

# get total number of visited locations
mtot_locs <- (mdf %>% 
               mutate(rlon = round(lon, -2), 
                      rlat = round(lat, -2)) %>% 
               dplyr::select(rlon, rlat) %>% 
               distinct() %>% 
               dim())[1]

# number of unique locations per species
ms <- mdf %>% 
  mutate(rlon = round(lon, -2), 
         rlat = round(lat, -2)) %>% 
  dplyr::select(rlon, rlat, sp_n) %>% 
  distinct(.keep_all = F)

# calculate prevalence for each species
m_sum <- ms %>% 
  group_by(sp_n) %>% 
  summarise(un_locs = length(rlon),
            prev = un_locs/mtot_locs) %>% 
  ungroup()
m_sum

ggplot(m_sum, aes(x = reorder(sp_n, -prev), y = prev)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle('Moth')

