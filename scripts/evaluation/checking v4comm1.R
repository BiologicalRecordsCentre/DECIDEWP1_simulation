

library(viridis)
library(tidyverse)
library(randomForest)
library(dismo)
library(tidyverse)
library(purrr)
library(cowplot)


# check predictions
sp50 <- read.csv('outputs/Sp50_asv1_v4community_1_50_sim_model_averages.csv')

sp50 %>% 
  # filter(method == 'initial' | method == 'none') %>% 
  ggplot(aes(x,y,fill = mean)) +
  geom_tile() +
  facet_wrap(~method) +
  scale_fill_viridis_c()


sp32 <- read.csv('outputs/Sp32_asv1_v4community_1_50_sim_model_averages.csv')


sp32 %>% 
  # filter(method == 'initial' | method == 'none') %>% 
  ggplot(aes(x,y,fill = mean)) +
  geom_tile() +
  facet_wrap(~method) +
  scale_fill_viridis_c()


# check results
res <- read.csv('outputs/asv1_v4community_1_50_sim_evaluation_table_alt.csv')
head(res)

res %>% 
  ggplot(aes(x = method, y = auc)) +
  geom_boxplot()

res %>% 
  mutate(d_auc = init_auc-auc,
         d_mse = init_mse-mse, 
         d_cor = init_corr-corr) %>% 
  pivot_longer(cols = c(d_auc, d_mse, d_cor)) %>% 
  ggplot() +
  geom_boxplot(aes(x = method, y = value, fill = name)) +
  ylab('initial - AS method')


res %>% 
  mutate(d_auc = init_auc-auc,
         d_mse = init_mse-mse, 
         d_cor = init_corr-corr) %>% 
  pivot_longer(cols = c(d_auc, d_mse, d_cor)) %>% 
  filter(method == 'none') %>% 
  ggplot() +
  geom_boxplot(aes(x = method, y = value, fill = name)) +
  ylab('initial - AS method')



## try looking at individual models
# are all models using same variables?

load('outputs/comm4_asv1_investigating/asv1_v4gam_SDMs_GBnew_Sp50_initial_AS_coverage.rdata')
(model_output$sdm_output[[1]])

load('outputs/comm4_asv1_investigating/asv1_v4lr_SDMs_GBnew_Sp50_initial_AS_coverage.rdata')
(model_output$sdm_output[[1]])


load('outputs/comm4_asv1_investigating/asv1_v4lr_SDMs_GBnew_Sp5_initial_AS_coverage.rdata')
(model_output$sdm_output[[1]])

load('outputs/comm4_asv1_investigating/asv1_v4gam_SDMs_GBnew_Sp5_initial_AS_coverage.rdata')
(model_output$sdm_output[[1]])


load('outputs/comm4_asv1_investigating/asv1_v4lr_SDMs_GBnew_Sp24_initial_AS_coverage.rdata')
(model_output$sdm_output[[1]])

load('outputs/comm4_asv1_investigating/asv1_v4gam_SDMs_GBnew_Sp24_initial_AS_coverage.rdata')
(model_output$sdm_output[[1]])

## Not entirely, GAMs seem to drop one variable sometimes why? Maybe variable selection?
#### Because of the dropping variables with too few unique data points for the number of knots!!!!


## What about new locations - do they make sense?
method = c("none", "uncertainty", "prevalence", "unc_plus_prev", "unc_plus_recs", "coverage")

method_locs <- data.frame()

for(m in seq_along(method)) {
  print(m)
  
  cov <- readRDS(paste0('outputs/comm4_asv1_investigating/asv2_v4community_1_50_sim_initial_AS_', method[m], '.rds'))
  
  obvs <- data.frame()
  
  for(i in 1:length(cov)) {
    obvs <- rbind(obvs, cbind(cov[[i]]$observations[,1:2], species = paste0('Sp',i))) 
  }
  
  obvs$id <- paste(obvs$lon, obvs$lat)
  obvs$method <- method[m]
  
  method_locs <- rbind(method_locs,obvs)
  
}

# get initial observations for each species
obvsinit <- readRDS('outputs/comm4_asv1_investigating/v4community_1_50_sim_initial.rds')

init_obvs <- data.frame()
for(i in 1:length(obvsinit)) {
  init_obvs <- rbind(init_obvs, cbind(obvsinit[[i]]$observations[,1:2], species = paste0('Sp',i))) 
}

obvsinit2 <- init_obvs %>% 
  dplyr::select(!species) %>% 
  distinct() %>% 
  mutate(method = 'initial',
         id = paste(lon, lat))
head(obvsinit2)


# find the new locs that weren't the initial locs
mlocs <- method_locs[!method_locs$id %in% obvsinit2$id,] %>% 
  dplyr::select(!species) %>% 
  distinct()
head(mlocs)

# check that there are no more than 2000 locations per method 
mlocs %>% group_by(method) %>% tally

# comine original and new locs
mlocs_all <- rbind(mlocs, obvsinit2)

ggplot() +
  geom_tile(data = sp32, aes(x,y,fill = mean)) +
  geom_point(data = mlocs_all, aes(lon,lat), size = 0.4, colour = 'red') +
  facet_wrap(~method) +
  scale_fill_viridis_c() +
  theme_bw()



### exploring cell weights
fls <- list.files('outputs/comm4_asv1_investigating/', 
                  pattern = 'cellweights', full.names = T)

fls_names <- list.files('outputs/comm4_asv1_investigating/', 
                        pattern = 'cellweights')

cell_weights <- data.frame()

for(x in 1:length(fls)) {
  
  cw <- read.csv(fls[x])
  
  cw$method <- gsub("_asv1_v4community_1_50_sim_initial_cellweights.csv",
                    '',
                    fls_names[x])
  
  if(any(cw$method=='coverage')) cw$cell_weights <- cw$butterfly_1km_effort
  
  cw2 <- cw %>% 
    dplyr::select(x, y, cell_weights, method)
  
  cell_weights <- rbind(cell_weights, cw2)
  
}


# ggplot() + 
#   geom_tile(data=cell_weights, aes(x,y,fill = log(cell_weights))) +
#   # geom_point(data = subset(mlocs_all, method == 'unc_plus_recs'),aes(lon,lat)) +
#   scale_fill_viridis_c() +
#   facet_wrap(~method)

head(cell_weights)

p <- cell_weights %>% 
  group_split(method) %>% 
  map(
    ~ggplot(., aes(x, y, fill = log(cell_weights))) + 
      geom_tile() + 
      scale_fill_viridis_c()+
      theme_bw() +
      theme(legend.position = 'none') +
      # xlim(0,680000) +
      facet_grid(~ method, labeller = function(x) label_value(x, multi_line = FALSE))
  ) %>% 
  plot_grid(plotlist = ., align = 'hv', axis='rl', ncol = 3)
p

# ggsave(p, 
#        filename = 'outputs/comm4_asv1_investigating/plots/v4comm1_cellweights_spatial_distrib.jpeg',
#        width = 9, height = 10)

head(cell_weights)
head(mlocs_all)

mlocs_all2 <- mlocs_all %>% 
  rename(x = lon, 
         y = lat, 
         cell_weights = id) %>% 
  mutate(cell_weights=0,
         type = 'obvs') %>% 
  filter(method != 'initial')

df <- rbind(cell_weights %>% mutate(type='cells'), mlocs_all2)

p2 <- df %>% 
  group_split(method) %>% 
  map(
    ~ggplot() + 
      geom_tile(data=subset(., type == 'cells'), aes(x, y, fill = log(cell_weights))) + 
      geom_point(data=subset(., type == 'obvs'), aes(x, y), colour = 'red', size = 0.7) + 
      scale_fill_viridis_c() +
      theme_bw() +
      theme(legend.position = 'none') +
      # xlim(0,680000) +
      facet_grid(~ method, labeller = function(x) label_value(x, multi_line = FALSE))
  ) %>% 
  plot_grid(plotlist = ., align = 'hv', axis='rl', ncol = 3)
p2


# ggsave(p2, 
#        filename = 'outputs/comm4_asv1_investigating/plots/v4comm1_cellweights_newobvs_spatial_distrib.jpeg',
#        width = 9, height = 10)
# 
# 
# ggsave(p2, 
#        filename = 'outputs/plots/v4comm1_cellweights_newobvs_spatial_distrib.jpeg',
#        width = 11, height = 10)




#################################################
###  Compare different AS versions (uptakes)  ###
#################################################

## community 2 as well as 1!

## What about new locations - do they make sense?
method = c("none", "uncertainty", "prevalence", "unc_plus_prev", "unc_plus_recs", "coverage")

method_locs <- data.frame()

for(com in 1:2) { # community
  
  for(a in 1:4){ # adaptive sampling version
    
    for(m in seq_along(method)) {
      print(m)
      
      cov <- readRDS(paste0('outputs/comm4_asv1_investigating/asv',a,'_v4community_', com, '_50_sim_initial_AS_', method[m], '.rds'))
      
      obvs <- data.frame()
      
      for(i in 1:length(cov)) {
        obvs <- rbind(obvs, cbind(cov[[i]]$observations[,1:2], species = paste0('Sp',i))) 
      }
      
      obvs$id <- paste(obvs$lon, obvs$lat)
      obvs$method <- method[m]
      obvs$asv <- a
      obvs$community <- com
      
      method_locs <- rbind(method_locs,obvs)
      
    }
  }
  
}

head(method_locs)


# get initial observations for each species

## community 1
obvsinit <- readRDS('outputs/comm4_asv1_investigating/v4community_1_50_sim_initial.rds')

init_obvs <- data.frame()
for(i in 1:length(obvsinit)) {
  init_obvs <- rbind(init_obvs, cbind(obvsinit[[i]]$observations[,1:2], species = paste0('Sp',i))) 
}

obvsinit <- init_obvs %>% 
  dplyr::select(!species) %>% 
  distinct() %>% 
  mutate(id = paste(lon, lat),
         method = 'initial',
         asv = 'init', 
         community = 1)
head(obvsinit)

meth_locsC1 <- subset(method_locs, community == 1)

# find the new locs that weren't the initial locs
mlocsC1 <- meth_locsC1[!meth_locsC1$id %in% obvsinit$id,] %>% 
  dplyr::select(!species) %>% 
  distinct()
head(mlocsC1)

## community 2
obvsinit2 <- readRDS('outputs/comm4_asv1_investigating/v4community_2_50_sim_initial.rds')

init_obvs2 <- data.frame()
for(i in 1:length(obvsinit2)) {
  init_obvs2 <- rbind(init_obvs2, cbind(obvsinit2[[i]]$observations[,1:2], species = paste0('Sp',i))) 
}

obvsinit2 <- init_obvs2 %>% 
  dplyr::select(!species) %>% 
  distinct() %>% 
  mutate(id = paste(lon, lat),
         method = 'initial',
         asv = 'init',
         community = 2)
head(obvsinit2)

meth_locsC2 <- subset(method_locs, community == 2)

# find the new locs that weren't the initial locs
mlocs2 <- meth_locsC2[!meth_locsC2$id %in% obvsinit2$id,] %>% 
  dplyr::select(!species) %>% 
  distinct()
head(mlocs2)

# comine original and new locs
mlocs_all <- rbind(mlocsC1, mlocs2, obvsinit, obvsinit2)

head(mlocs_all)

ggplot(subset(mlocs_all, asv!='init' & community == 1), aes(x=lon, y=lat)) +
  geom_point(size = 0.7) +
  facet_wrap(asv~method, ncol = 6)


ggplot(subset(mlocs_all, asv!='init'), aes(x=lon, y=lat)) +
  geom_point(size = 0.7) +
  facet_wrap(community~method, ncol = 6)




