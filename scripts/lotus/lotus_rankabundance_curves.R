
###   Get rank abundance curves for all communities

##########
##  This is run in an R terminal on lotus

# package for sorting
library(tidyverse)

# make sure we're in the right directory
setwd('/gws/nopw/j04/ceh_generic/thoval/DECIDE/simulation')

# get the directories
dirs <- config::get("LOTUSpaths")

# functions to use
source(paste0(dirs$inpath, 'reformat_simulated_data.R'))

# number of communities to get through
n_communities = 1:20

# lapply through all communities
all_comm <- lapply(n_communities, FUN = function(com_n){
  
  # read in community data
  spdata <- paste0(dirs$commpath, sprintf("community_%i_50_sim/community_%i_50_sim_initial.rds", com_n, com_n))
  
  # load in community data
  community <- readRDS(spdata)
  
  # extract the prevalences
  sp_prev <- as.data.frame(do.call(rbind, lapply(1:length(community), FUN=function(x) cbind(community_name = paste0('community', com_n), sp_name = paste0('sp', x), prev = community[[x]]$prevalence))))
  
  return(sp_prev)
  
})

# bind all together
all_comm_out <- do.call(rbind, all_comm)

# write to csv
write.csv(all_comm_out, file = paste0(dirs$outpath, 'combined_community_', paste0(range(n_communities), collapse = '_'), '_output.csv'))


################################### Home PC
# read csv on base computer
ra_c <- read.csv('Outputs/combined_community_1_20_output.csv')

ra_c <- ra_c %>% 
  group_by(community_name) %>% 
  arrange(-prev) %>% 
  mutate(rank = 1:50) %>% 
  arrange(X)

head(ra_c)

ggplot(ra_c, aes(x = rank, y = prev, colour = community_name)) +
  geom_line() +
  ylim(0,0.5) +
  theme_classic() +
  ylab('Prevalence')
