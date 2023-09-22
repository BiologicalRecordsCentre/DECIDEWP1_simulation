
## Supplementary material figures

library(tidyverse)

write = TRUE


# Butterfly rank abundance curve
{
  bdf <- read_csv("../DECIDE_WP1/data/edited_insect_data/butterfly/butterfly_EastNorths_no_duplicates_2021_12_06.csv") %>% 
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
  
  rankabund_butterfly <- ggplot(b_sum, aes(x = reorder(sp_n, -prev), y = prev)) +
    geom_point() +
    theme_classic() +
    ylab('Prevalence') +
    xlab('') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12)) +
    ggtitle('Butterfly')
  rankabund_butterfly
  
  if(write){
    ggsave(rankabund_butterfly, filename = 'outputs/plots/paper/supp_matt_rankabund_butterfly.png',
           width = 9, height = 5)
  }  
}

# communities rank abundance curves
{
  ra_c <- read.csv('outputs/plots/paper/v4combined_community_1_50_output.csv')
  
  ra_c <- ra_c %>% 
    group_by(community_name) %>% 
    arrange(-prev) %>% 
    mutate(rank = 1:50) %>% 
    arrange(X)
  
  head(ra_c)
  
  rank_abund_comm <- ggplot(ra_c, aes(x = rank, y = prev, colour = community_name)) +
    geom_line() +
    ylim(0,0.5) +
    theme_classic() +
    ylab('Prevalence') +
    xlab('Rank') + 
    labs(colour = 'Community name')
  rank_abund_comm
  
  if(write){
    ggsave(rank_abund_comm, filename = 'outputs/plots/paper/supp_matt_rankabund_communities_all.png',
           width = 9, height = 5)
  }  
}