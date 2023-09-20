

library(tidyverse)
library(viridis)
library(ggridges)
library(patchwork)

write = TRUE


## set percentage change
p_c <- 5


## sort data out
{
  
  # for changing labels
  meth_names <- list(
    "initial",
    "Business\nas usual",
    "Gap-filling",
    "Rare species",
    "Uncertainty only",
    "Uncertainty of\nrare species",
    "Gap-filling\nwith uncertainty"
  )
  
  ## load uptake values
  obvs_df <- read.csv('outputs/v4Community/v4n_new_obvs_perc_increase_1_50_spp50.csv', 
                      stringsAsFactors = FALSE)
  
  # load each of the evaluation files 
  cdf_0.1_uptake <- read.csv('outputs/v4Community/asv1_v4combined_outputs_comm1_50_spp50_v2.csv', 
                             stringsAsFactors = FALSE) %>%
    mutate(uptake = '0.1',
           asv = 'asv1')
  cdf_0.01_uptake <- read.csv('outputs/v4Community/asv2_v4combined_outputs_comm1_50_spp50_v2.csv', 
                              stringsAsFactors = FALSE) %>%
    mutate(uptake = '0.01',
           asv = 'asv2')
  
  cdf_0.5_uptake <- read.csv('outputs/v4Community/asv4_v4combined_outputs_comm1_50_spp50_v2.csv', 
                             stringsAsFactors = FALSE) %>%
    mutate(uptake = '0.5',
           asv = 'asv4')
  
  cdf_0_uptake <- read.csv('outputs/v4Community/asv3_v4combined_outputs_comm1_50_spp50_v2.csv', 
                           stringsAsFactors = FALSE) %>%
    mutate(uptake = '0',
           asv = 'asv3')
  
  # combine all the files
  cdf <- rbind(cdf_0_uptake, cdf_0.1_uptake,cdf_0.01_uptake,cdf_0.5_uptake)
  
  # calculate differences
  init_tab <- cdf[cdf$method =='initial',]
  colnames(init_tab) <- paste0('initial_', colnames(init_tab))
  et <- cdf
  
  et$init_mse <- init_tab$initial_mse[match(paste0(et$species, et$community),
                                            paste0(init_tab$initial_species, init_tab$initial_community))]
  
  et$init_medse <- init_tab$initial_medianse[match(paste0(et$species, et$community),
                                                   paste0(init_tab$initial_species, init_tab$initial_community))]
  et$init_corr <- init_tab$initial_corr[match(paste0(et$species, et$community),
                                              paste0(init_tab$initial_species, init_tab$initial_community))]
  et$init_auc <- init_tab$initial_auc[match(paste0(et$species, et$community),
                                            paste0(init_tab$initial_species, init_tab$initial_community))]
  et$init_mean_sd <- init_tab$initial_mean_sd[match(paste0(et$species, et$community),
                                                    paste0(init_tab$initial_species, init_tab$initial_community))]
  et$init_max_sd <- init_tab$initial_max_sd[match(paste0(et$species, et$community),
                                                  paste0(init_tab$initial_species, init_tab$initial_community))]
  
  
  ## get differences
  et <- et %>% 
    mutate(delta_mse = mse - init_mse,
           delta_medse = medianse - init_medse,
           delta_corr = corr - init_corr,
           delta_auc = auc - init_auc,
           delta_mean_sd = mean_sd - init_mean_sd,
           delta_max_sd = max_sd - init_max_sd)
  
  ## average across communities
  comm_df <- et %>%
    group_by(community, method, uptake) %>%
    summarise(mse = mean(mse, na.rm = TRUE),
              medse = mean(medianse, na.rm = TRUE),
              corr = mean(corr, na.rm = TRUE),
              auc = mean(auc, na.rm = TRUE),
              mean_sd = mean(mean_sd, na.rm = TRUE),
              max_sd = max(max_sd, na.rm = TRUE),
              prev = median(prevalence, na.rm = TRUE),
              init_mse = mean(init_mse, na.rm = TRUE),
              init_medse = mean(init_medse, na.rm = TRUE),
              init_corr = mean(init_corr, na.rm = TRUE),
              init_auc = mean(init_auc, na.rm = TRUE),
              init_mean_sd = mean(init_mean_sd, na.rm = TRUE),
              init_max_sd = max(init_max_sd, na.rm = TRUE),
              prev = median(prev)) %>%
    ungroup() %>% 
    mutate(delta_mse = sqrt(mse) - sqrt(init_mse),
           delta_medse = (medse) - (init_medse),
           delta_corr = corr - init_corr,
           delta_auc = auc - init_auc,
           delta_mean_sd = mean_sd - init_mean_sd,
           delta_max_sd = max_sd - init_max_sd,
           method = factor(method, 
                           levels=c("initial", "none", "coverage", "prevalence",  
                                    "uncertainty", "unc_plus_prev", "unc_plus_recs")))
  
  
  ## load all the observations for a given AS version
  asv1 <- read_csv('outputs/v4Community/asv1_v4all_observations.csv') %>% 
    mutate(uptake = 0.1,
           asv = 'asv1')
  asv2 <- read_csv('outputs/v4Community/asv2_v4all_observations.csv') %>% 
    mutate(uptake = 0.01,
           asv = 'asv2')
  asv3 <- read_csv('outputs/v4Community/asv3_v4all_observations.csv') %>% 
    mutate(uptake = 0,
           asv = 'asv3')
  asv4 <- read_csv('outputs/v4Community/asv4_v4all_observations.csv') %>% 
    mutate(uptake = 0.5,
           asv = 'asv4')
  
  df <- rbind(asv1,asv2,asv3,asv4)
  
  # initial locations
  init_nona <- df[which(df$method == 'initial'),]
  
  # all others
  meths <- df#[df$method != 'initial',]
  
  # remove initial locations - only looking at new locations
  new_locs <- meths[!meths$id %in% init_nona$id,] ## id is lon lat id
  
  # initial locations for each species
  sp_init_locs <- init_nona %>% 
    # na.omit() %>% 
    group_by(method, community, species, prevalence, asv, uptake) %>% 
    summarise(n = sum(Observed, na.rm = TRUE)) %>% 
    mutate(id = paste(community, species, prevalence, asv, uptake, sep = '_'))
  
  # new locations for each species
  sp_new_locs <- new_locs %>% 
    # na.omit() %>% 
    group_by(method, community, species, prevalence, asv, uptake) %>% 
    summarise(n = sum(Observed)) %>% 
    mutate(id = paste(community, species, prevalence, asv, uptake, sep = '_'))
  
  # bind the two together
  init_new_locs <- sp_new_locs %>% 
    rowwise() %>% 
    mutate(n_init_obs = sp_init_locs$n[match(id, sp_init_locs$id)],
           id2 = paste(community, method, species, prevalence, asv, uptake, sep = '_'))
  
  ### get model improvements above a certain percentage
  etp <- et %>% 
    # group_by(method, uptake) %>% 
    mutate(prev_cat = dplyr::ntile(prevalence, 10),
           auc_cat = dplyr::ntile(init_auc, 10),
           mse_cat = dplyr::ntile(init_mse, 10),
           medse_cat = dplyr::ntile(init_medse, 10),
           corr_cat = dplyr::ntile(init_corr, 10)) %>% 
    # ungroup() %>% 
    # rowwise() %>% 
    mutate(perc_inc_auc = (delta_auc)/(init_auc)*100,
           perc_inc_corr = (delta_corr)/(init_corr)*100,
           perc_inc_mse = (delta_mse)/(init_mse)*100,
           perc_inc_medse = (delta_medse)/(init_medse)*100,
           perc_imp_auc = ifelse(perc_inc_auc>= p_c, p_c, 
                                 ifelse(perc_inc_auc<= -p_c, -p_c, 0)),
           perc_imp_corr = ifelse(perc_inc_corr>=p_c, p_c, 
                                  ifelse(perc_inc_corr<= -p_c, -p_c, 0)),
           perc_imp_mse = ifelse(perc_inc_mse<= -p_c, p_c, 
                                 ifelse(perc_inc_mse>= p_c, -p_c, 0)),
           perc_imp_medse = ifelse(perc_inc_medse<= -p_c, p_c, 
                                   ifelse(perc_inc_medse>= p_c, -p_c, 0)),
           method = factor(method, 
                           levels=c("initial", "none", "coverage", "prevalence", 
                                    "uncertainty", "unc_plus_prev", "unc_plus_recs")))
  
  # number of models with > x% increase in each community
  nmods <- etp %>%
    group_by(community, method, uptake) %>%
    summarise(n_mods_auc_1 = sum(perc_inc_auc>1, na.rm = TRUE),
              n_mods_auc_2 = sum(perc_inc_auc>2, na.rm = TRUE),
              n_mods_auc_5 = sum(perc_inc_auc>5, na.rm = TRUE),
              n_mods_mse_1 = sum(perc_inc_mse< -1, na.rm = TRUE),
              n_mods_mse_2 = sum(perc_inc_mse< -2, na.rm = TRUE),
              n_mods_mse_5 = sum(perc_inc_mse< -5, na.rm = TRUE),
              n_mods_medse_1 = sum(perc_inc_medse< -1, na.rm = TRUE),
              n_mods_medse_2 = sum(perc_inc_medse< -2, na.rm = TRUE),
              n_mods_medse_5 = sum(perc_inc_medse< -5, na.rm = TRUE),
              n_mods_corr_1 = sum(perc_inc_corr>1, na.rm = TRUE),
              n_mods_corr_2 = sum(perc_inc_corr>2, na.rm = TRUE),
              n_mods_corr_5 = sum(perc_inc_corr>5, na.rm = TRUE)) 
  
  nmods_l <- pivot_longer(nmods, cols = 4:15) %>% 
    rowwise() %>% 
    mutate(eval_type = ifelse(grepl(x = name, pattern = 'auc'), 'auc',
                              ifelse(grepl(x = name, pattern = 'mse'), 'mse', 
                                     ifelse(grepl(x = name, pattern = 'corr'), 'corr',
                                            ifelse(grepl(x = name, pattern = 'medse'), 'medse', 'WRONG')))),
           inc_amount = as.numeric(gsub("[^\\d]+", "", name, perl=TRUE)),
           method = factor(method, 
                           levels=c("initial", "none", "coverage", "prevalence", 
                                    "uncertainty", "unc_plus_prev", "unc_plus_recs")))
  
  
  # data for alternate figure 4
  
  etp2 <- et %>% 
    # group_by(method, uptake) %>% 
    mutate(prev_cat = dplyr::ntile(prevalence, 10),
           auc_cat = dplyr::ntile(init_auc, 10),
           mse_cat = dplyr::ntile(init_mse, 10),
           medse_cat = dplyr::ntile(init_medse, 10),
           corr_cat = dplyr::ntile(init_corr, 10)) %>% 
    # ungroup() %>% 
    # rowwise() %>% 
    mutate(perc_inc_auc = (delta_auc)/(init_auc)*100,
           perc_inc_corr = (delta_corr)/(init_corr)*100,
           perc_inc_mse = (delta_mse)/(init_mse)*100,
           perc_inc_medse = (delta_medse)/(init_medse)*100,
           perc_imp_auc = ifelse(perc_inc_auc>= p_c, p_c, 
                                 ifelse(perc_inc_auc<= -p_c, -p_c, 0)),
           perc_imp_corr = ifelse(perc_inc_corr>=p_c, p_c, 
                                  ifelse(perc_inc_corr<= -p_c, -p_c, 0)),
           perc_imp_mse = ifelse(perc_inc_mse<= -p_c, p_c, 
                                 ifelse(perc_inc_mse>= p_c, -p_c, 
                                        ifelse(perc_inc_mse>= -p_c & perc_inc_mse< 0, 2.5,
                                               ifelse(perc_inc_mse<= p_c & perc_inc_mse> 0, -2.5, 0)))),
           perc_imp_medse = ifelse(perc_inc_medse<= -p_c, p_c, 
                                 ifelse(perc_inc_medse>= p_c, -p_c, 
                                        ifelse(perc_inc_medse> -p_c & perc_inc_medse< 0, 2.5,
                                               ifelse(perc_inc_medse< p_c & perc_inc_medse> 0, -2.5, 0)))),
           method = factor(method, 
                           levels=c("initial", "none", "coverage", "prevalence",  
                                    "uncertainty", "unc_plus_prev", "unc_plus_recs")))
  
  
  etp_p2 <- subset(etp2, method != 'initial')
  levels(etp_p2$method) <- unlist(meth_names)
  
  # load the probability sampling layers
  {
    
    source('scripts/functions/filter_distance.R')
    
    ### exploring cell weights
    fls <- list.files('outputs/comm4_asv1_investigating/', 
                      pattern = 'cellweights', full.names = TRUE)
    
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
    
    method_unique <- unique(cell_weights$method)
    
    cwt_rast <- list()
    for(md in seq_along(method_unique)) {
      
      cell_wt_method <- cell_weights[cell_weights$method == method_unique[md],1:3]
      cwt_rast[[md]] <- raster::rasterFromXYZ(cell_wt_method)
      
    }
    
    cwt_rast[[1]] <- raster::crop(cwt_rast[[1]], cwt_rast[[2]])  
    
    all_comb <- raster::stack(cwt_rast)
    names(all_comb) <- method_unique  
    
    crp_rst <- filter_distance(all_comb,
                               location = c(-2.2, 53.8),
                               distance = 20000)
    
    as_layers <- as.data.frame(crp_rst, xy = TRUE) %>% 
      pivot_longer(cols = 3:8)
    
  }
  
  ## load the observations 
  {
    
    #################################################
    ###  Compare different AS versions (uptakes)  ###
    #################################################
    
    ## community 2 as well as 1!
    
    ## What about new locations - do they make sense?
    method = c("none", "uncertainty", "prevalence", "unc_plus_prev", "unc_plus_recs", "coverage")
    
    method_locs <- data.frame()
    
    for(com in 1) { # community
      
      for(a in 1){ # adaptive sampling version
        
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
    
    # comine original and new locs
    mlocs_all <- rbind(mlocsC1, obvsinit)
    
    mlocs_all_sub <- mlocs_all[!is.na(raster::extract(crp_rst[[2]], mlocs_all[,c(1,2)])), ]
    
  }
  
  # environmental data
  env_data <- raster::stack('data/environmental_data/envdata_1km_no_corr_noNA.grd')[[22:24]]
  
  library(sf)
  library(stars)
  library(rnaturalearth)
  
  # get world map
  world <- ne_countries(scale = "medium", type = 'map_units', returnclass = "sf")
  class(world)
  
  # get czech map
  uk <- world[world$admin == 'United Kingdom',]
  gb <- uk[uk$name != 'N. Ireland',] %>% 
    st_transform(crs = 27700)
  
  # get community data
  c4 <- readRDS('outputs/v4Community/v4community_1_50_sim_initial.rds')
  c4cov <- readRDS('outputs/comm4_asv1_investigating/asv1_v4community_1_50_sim_initial_AS_uncertainty.rds')
  load('outputs/comm4_asv1_investigating/asv1_v4lr_SDMs_GBnew_Sp5_initial_AS_uncertainty.rdata')
  
  
  
  ## determine differences in observations/characteristics per community/uptake
  loc_summ <- new_locs %>% 
    na.omit() %>% 
    mutate(locid = paste(lon, lat, sep = '_')) %>% 
    group_by(method, community, uptake) %>% # work out some community-level summaries
    summarise(total_obs = sum(Observed),
              unique_locs = length(unique(locid)),
              av_obs_per_loc = sum(Observed)/unique_locs, # number of observations per location
              diversity = length(unique(species)), # n unique species
              av_div_per_loc = sum(length(unique(species)))/unique_locs) %>% # n unique species per location
    mutate(method = factor(method, 
                           levels=c("initial", "none", "coverage", "prevalence",  
                                    "uncertainty", "unc_plus_prev", "unc_plus_recs")))
  
  prev_per_loc <- new_locs %>% 
    na.omit() %>% 
    mutate(locid = paste(lon, lat, sep = '_')) %>% 
    group_by(method, community, uptake, locid) %>% # work out median prevalence per location
    summarise(prev_per_loc = median(prevalence)) %>% 
    mutate(method = factor(method, 
                           levels=c("initial", "none", "coverage", "prevalence",  
                                    "uncertainty", "unc_plus_prev", "unc_plus_recs")))
  
  
  # calculate prevalence for each community/method/uptake
  prev_df <- new_locs[, c('method', 'community', 'uptake', 'species', 'prevalence')]
  prev_df_unique <- prev_df[!duplicated(prev_df),]
  prev_df_unique <- prev_df_unique %>% 
    na.omit() %>% 
    group_by(method, community, uptake) %>% 
    summarise(med_prev = median(prevalence)) %>% 
    mutate(method = factor(method, 
                           levels=c("initial", "none", "coverage", "prevalence",  
                                    "uncertainty", "unc_plus_prev", "unc_plus_recs")))
  
  # rename methods
  levels(loc_summ$method) <- unlist(meth_names)
  levels(prev_df_unique$method) <- unlist(meth_names)
  levels(prev_per_loc$method) <- unlist(meth_names)
  
  # # write out file to create 'recorder preferences' plots
  # write.csv(loc_summ, file = 'outputs/plots/paper/recprefs_loc_summ.csv')
  
  
}



## plotting -----

## figure 2 - plots for the AS methodological figure
head(as_layers)
as_method <- method[2:5]

for(md in 1:length(as_method)) {
  
  print(as_method[md])
  
  p <- ggplot(subset(as_layers, name == as_method[md]), 
              aes(x,y,fill = log(value))) +
    geom_raster()+
    scale_fill_viridis(na.value="transparent") +
    theme_void() +
    theme(legend.position = 'none') 
  print(p)
  if(write){
    ggsave(p, filename = paste0('outputs/plots/paper/figure2_ASmethodplot_', as_method[md], '.png'),
           bg = 'transparent', width = 2.2, height = 2.2)
  }
}

for(m in method[c(1,6)]){
  p_nas <- ggplot() +
    geom_raster(data = subset(as_layers, name == 'coverage'), aes(x,y,fill = log(value))) +
    geom_point(data = subset(mlocs_all_sub, method == 'initial'), aes(lon,lat), colour = 'black') +
    geom_point(data = subset(mlocs_all_sub, method == m), aes(lon,lat), colour = 'red') +
    # scale_colour_manual(values = c('black', 'red')) +
    scale_fill_viridis(na.value = 'transparent') +
    theme_void() +
    theme(legend.position = 'none')
  # print(p_nas)
  
  if(write){
    ggsave(p_nas, filename = paste0('outputs/plots/paper/figure2_ASmethodplot_', m, '.png'),
           bg = 'transparent', width = 2.2, height = 2.2)
  }
}


## figure 1 - plots for the AS methodological figure
## plot of observation data
# sample to mimic only a few species
uncert <- all_comb[['uncertainty']] %>% 
  as.data.frame(xy = TRUE) %>% 
  mutate(uncert = ifelse(is.na(uncertainty), NA, 1))

initobs <- subset(mlocs_all, asv == 'init')

out_obs <- data.frame()
for(n in c(200, 600, 1000)) {
  samp_obvs_ind <- sample(1:nrow(initobs), size = n)
  
  out_obs <- rbind(out_obs, cbind(initobs[samp_obvs_ind,], species = paste0('sp',n)))
}

p_obs <- ggplot() +
  geom_sf(data = gb, fill = 'white', col = 'white', size = 5) +
  scale_fill_manual(na.value = 'transparent') +
  geom_point(data = out_obs, aes(lon, lat, colour = species), size = 1) +
  scale_colour_manual(values = c("#E69F00", "#009E73", "#D55E00")) +
  theme_void() +
  theme(legend.position = 'none')
p_obs

if(write){
  ggsave(p_obs, filename = paste0('outputs/plots/paper/figure1_simmethodplot_sp_obs.png'),
         bg = 'transparent', width = 4.46, height = 7)
}

## environmental data
edf <- as.data.frame(env_data, xy = TRUE)
head(edf)


ed1 <- ggplot() +
  geom_sf(data = gb, fill = 'white', col = 'white', size = 4) +
  geom_tile(data = edf, aes(x,y, fill = AnnualTemp)) +
  scale_fill_continuous(na.value = 'transparent') +
  theme_void() +
  theme(legend.position = 'none')

ed2 <- ggplot() +
  geom_sf(data = gb, fill = 'white', col = 'white', size = 4) +
  geom_tile(data = edf, aes(x,y, fill = MeanDiRange)) +
  scale_fill_viridis(na.value = 'transparent') +
  theme_void() +
  theme(legend.position = 'none')

ed3 <- ggplot() +
  geom_sf(data = gb, fill = 'white', col = 'white', size = 4) +
  geom_tile(data = edf, aes(x,y, fill = Isotherm)) +
  scale_fill_viridis(na.value = 'transparent', option = 'C') +
  theme_void() +
  theme(legend.position = 'none')

if(write) {  
  ggsave(ed1, filename = paste0('outputs/plots/paper/figure1_simmethodplot_antemp.png'),
         bg = 'transparent', width = 4.46, height = 7)  
  ggsave(ed2, filename = paste0('outputs/plots/paper/figure1_simmethodplot_mnditemp.png'),
         bg = 'transparent', width = 4.46, height = 7)  
  ggsave(ed3, filename = paste0('outputs/plots/paper/figure1_simmethodplot_isotherm.png'),
         bg = 'transparent', width = 4.46, height = 7)
  
}


## true prob of occurrence 
c4sp3 <- c4[[5]]$pres_abs %>% 
  as.data.frame(xy=TRUE)
c4covsp3 <- c4cov[[5]]$pres_abs %>% 
  as.data.frame(xy=TRUE)


true_sp <- ggplot() +
  geom_sf(data = gb, fill = 'white', col = 'white', size = 5) +
  geom_raster(data = c4sp3, aes(x,y,fill = layer)) +
  scale_fill_manual(na.value = 'transparent', values = c('grey', 'green')) +
  theme_void() +
  theme(legend.position = 'none')

pred_sp1 <- model_output$predictions
load('outputs/comm4_asv1_investigating/asv1_v4lr_SDMs_GBnew_Sp50_initial_AS_uncertainty.rdata')
pred_sp2 <- model_output$predictions


pred_sp_pl <- ggplot() +
  geom_sf(data = gb, fill = 'white', col = 'white', size = 4) +
  geom_raster(data = pred_sp1, aes(x,y,fill = mean)) +
  scale_fill_viridis(na.value = 'transparent') +
  theme_void() +
  theme(legend.position = 'none')
pred_sp_pl


pred_sp_pl2 <- ggplot() +
  geom_sf(data = gb, fill = 'white', col = 'white', size = 4) +
  geom_raster(data = pred_sp2, aes(x,y,fill = mean)) +
  scale_fill_viridis(na.value = 'transparent') +
  theme_void() +
  theme(legend.position = 'none')
pred_sp_pl2


if(write) {  
  ggsave(true_sp, filename = paste0('outputs/plots/paper/figure1_simmethodplot_truepres.png'),
         bg = 'transparent', width = 4.46, height = 7)  
  
  ggsave(pred_sp_pl, filename = paste0('outputs/plots/paper/figure1_simmethodplot_predpres_sp1.png'),
         bg = 'transparent', width = 4.46, height = 7)  
  
  ggsave(pred_sp_pl2, filename = paste0('outputs/plots/paper/figure1_simmethodplot_predpres_sp2.png'),
         bg = 'transparent', width = 4.46, height = 7)  
  
  
}

## figure 3 - model improvements + N models > %

sd(comm_df$delta_mse)  
hist(comm_df$delta_mse)


cno1 <- ggplot(comm_df[comm_df$method!='initial' & comm_df$uptake!=0,] %>% 
                 mutate(facet = 'MSE'), 
               aes(x=method, y = delta_mse, fill = factor(uptake))) +
  geom_boxplot() +
  theme_bw() + 
  ylim(-(2*sd(comm_df$delta_mse)), 2*sd(comm_df$delta_mse)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  xlab('') +
  ylab('Delta MSE (lower = better)') +
  scale_fill_manual(name = 'Uptake (%)', labels = c(1, 10, 50),
                    values = c("#E69F00", "#56B4E9", "#009E73")) +
  # facet_wrap(~facet) +
  scale_x_discrete(labels= c('Business \n as usual', 'Gap-filling', 'Rare species',
                             'Uncertainty only', 'Uncertainty of \n rare species', 
                             'Gap-filling \n with uncertainty')) +
  theme(text = element_text(size = 12),
        axis.text.x = element_blank()) #element_text(size = 12, angle = 0, vjust = 0))

cno1

ggplot(comm_df[comm_df$method!='initial' & comm_df$uptake!=0,] %>% 
         mutate(facet = 'MSE'), 
       aes(x=method, y = delta_medse, fill = factor(uptake))) +
  geom_boxplot() 

msen <- ggplot(subset(nmods_l, uptake != 0 & inc_amount == 1 & method != 'initial' & eval_type == 'mse'), 
               aes(x = method, y = value, fill = factor(uptake))) +
  geom_boxplot() +
  scale_fill_manual(name = "Uptake (%)", labels = c(1, 10, 50),
                    values = c("#E69F00", "#56B4E9", "#009E73")) +
  # facet_wrap(~eval_type,ncol = 3) +
  ylab('Number of models with\n>1% improvement') +
  xlab('') +
  # ggtitle('Number of models with 1% improvement for different uptake values') +
  theme_bw() +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  scale_x_discrete(labels= c('Business \n as usual', 'Gap-filling', 'Rare species',
                             'Uncertainty\nonly', 'Uncertainty of \n rare species', 
                             'Gap-filling \n with uncertainty')) + 
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, vjust = 0))
msen


fig3 <- cno1 / msen + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'a')
fig3

if(write){
  ggsave(fig3, filename = 'outputs/plots/paper/figure3_improve_uptake.png',
         width = 9.5, height = 9.5)
}




## figure 4 -- proportion of models with >x% changes


etp_p <- subset(etp, method != 'initial')
levels(etp_p$method) <- unlist(meth_names)


fig4 <- ggplot(na.omit(subset(etp_p, uptake == 0.1 & method != 'initial')), 
               aes(x = prev_cat, fill = factor(perc_imp_mse))) +
  geom_bar(position="fill") +
  ylab('Proportion of models') +
  xlab('Prevalence category') +
  # ylim(0,0.25) +
  facet_grid(~method) +
  scale_fill_manual(name = 'Change in MSE (%)', labels = c('< -5', '0', '> 5'),
                    values = c("#E69F00", "#56B4E9", "#009E73")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1),
        text = element_text(size = 12))
fig4

if(write){
  ggsave(fig4, filename = 'outputs/plots/paper/figure4_propmods_0_5.png',
         width = 7, height = 4)
}


## figure 4 alternate --- proportion of models with more groups of change



fig4_s <- ggplot(na.omit(subset(etp_p2, uptake == 0.5 & method != 'initial')), 
                 aes(x = prev_cat, fill = factor(perc_imp_mse))) +
  geom_bar(position="fill") +
  ylab('Proportion of models') +
  xlab('Prevalence category') +
  # ylim(0,0.25) +
  facet_grid(~method) +
  scale_fill_manual(name = 'Improvement in\nmodel MSE (%)', 
                    # labels = c('< -5', '-5 to -1', 
                    #            '-1 to 1', '1 to 5', '> 5'),
                    labels = c('< -5', '-5 to 0', 
                               '0 to 5', '> 5'),
                    values = c("#E69F00", "#56B4E9", "#009E73",
                               "#0072B2")) +#, "#D55E00")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1),
        text = element_text(size = 14))
fig4_s

# median
ggplot(na.omit(subset(etp_p2, uptake == 0.5 & method != 'initial')), 
       aes(x = prev_cat, fill = factor(perc_imp_medse))) +
  geom_bar(position="fill") +
  ylab('Proportion of models') +
  xlab('Prevalence category') +
  # ylim(0,0.25) +
  facet_grid(~method) +
  scale_fill_manual(name = 'Improvement in model\nMSE (%)',
                    # labels = c('< -5', '-5 to -1',
                    #            '-1 to 1', '1 to 5', '> 5'),
                    labels = c('< -5', '-5 to 0', '0',
                               '0 to 5', '> 5'),
                    values = c("#E69F00", "#56B4E9", "#009E73",
                               "#0072B2", "#D55E00")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1),
        text = element_text(size = 14))

if(write){
  ggsave(fig4_s, filename = 'outputs/plots/paper/figure4_propmods_0_1_5.png',
         width = 11, height = 5)
}




## figure 5 --- Exploring benefits of AS

## plots
# total number of observations per community
s1 <- ggplot(data = subset(loc_summ, uptake==0.5), 
             aes(x=method, y=total_obs)) +
  geom_boxplot() +
  ylab('Total new observations\nper community') +
  xlab('') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        text = element_text(size = 12))

s2 <- ggplot(data = subset(loc_summ, uptake==0.5), 
             aes(x=method, y=diversity)) +
  geom_boxplot() +
  ylab('Species diversity\nper community') +
  xlab('') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        text = element_text(size = 12))

s3 <- ggplot(data = subset(prev_df_unique, uptake==0.5), 
             aes(x=method, y=med_prev)) +
  geom_boxplot() +
  ylab('Median prevalence\nper community') +
  xlab('') +
  theme_classic() +
  theme(text = element_text(size = 12))

fig5 <- s1/s2/s3 +
  plot_annotation(tag_levels = 'a')
fig5

if(write){
  ggsave(fig5, filename = 'outputs/plots/paper/figure5_benefits_as.png',
         width = 7, height = 8)
}


### plots per location ---

# total observations per location
s1.5 <- ggplot(data = subset(loc_summ, uptake==0.5),
               aes(x=method, y=av_obs_per_loc)) +
  geom_boxplot() +
  ylab('Observations\nper visit') +
  xlab('') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        text = element_text(size = 12))


# diverstiy per location
s2.5 <- ggplot(data = subset(loc_summ, uptake==0.5),
               aes(x=method, y=av_div_per_loc)) +
  geom_boxplot() +
  ylab('Unique species\nper visit') +
  xlab('') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        text = element_text(size = 12))


# diverstiy per location
s3.5 <- ggplot(data = subset(prev_per_loc, uptake==0.5),
               aes(x=method, y=prev_per_loc)) +
  geom_boxplot() +
  ylab('Prevalence\nper visit') +
  xlab('') +
  theme_classic() +
  theme(text = element_text(size = 12))


s1.5 / s2.5 / s3.5


# looking at number of things by grouping by location - not average!!
tdf <- new_locs %>% 
  na.omit() %>% 
  mutate(locid = paste(lon, lat, sep = '_')) %>% 
  group_by(method, community, uptake, locid) %>% 
  summarise(visit_numbers = length(locid),
            obs_per_loc = sum(Observed), # number of observations per location
            diversity = length(unique(species)),
            prev_per_loc = median(prevalence)) %>% 
  mutate(method = factor(method, 
                         levels=c("initial", "none", "coverage", "prevalence",  
                                  "uncertainty", "unc_plus_prev", "unc_plus_recs")),
         same_obs_div = ifelse(obs_per_loc == diversity, 1, 0))


levels(tdf$method) <- unlist(meth_names)


# total observations per location
s1.2 <- ggplot(data = subset(tdf, uptake==0.5),
               aes(x=method, y=visit_numbers)) +
  geom_boxplot() +
  ylab('num observations\nper visit') +
  xlab('') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        text = element_text(size = 12))

# total observations per location
s1.5 <- ggplot(data = subset(tdf, uptake==0.5),
               aes(x=method, y=obs_per_loc)) +
  geom_boxplot() +
  ylab('Total observations\nper visit') +
  xlab('') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        text = element_text(size = 12))


# diverstiy per location
s2.5 <- ggplot(data = subset(tdf, uptake==0.5),
               aes(x=method, y=diversity)) +
  geom_boxplot() +
  ylab('Unique species\nper visit') +
  xlab('') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        text = element_text(size = 12))


# diverstiy per location
s3.5 <- ggplot(data = subset(tdf, uptake==0.5),
               aes(x=method, y=prev_per_loc)) +
  geom_boxplot() +
  ylab('Prevalence\nper visit') +
  xlab('') +
  theme_classic() +
  theme(text = element_text(size = 12))


s1.2 / s1.5 / s2.5 / s3.5


ggplot(subset(tdf, uptake==0.5), aes(visit_numbers, diversity)) +
  geom_point() + 
  xlab('number of visits to a location') +
  facet_wrap(~method)



# %>% ungroup() %>% 
#   group_by(method, community, uptake) %>% # work out some community-level summaries
#   mutate(total_obs = sum(Observed),
#             unique_locs = length(unique(locid)),
#             av_obs_per_loc = sum(Observed)/unique_locs, # number of observations per location
#             diversity = length(unique(species)), # n unique species
#             av_div_per_loc = length(unique(species))/unique_locs) %>% # n unique species per location
#   mutate(method = factor(method, 
#                          levels=c("initial", "none", "coverage", "prevalence",  
#                                   "uncertainty", "unc_plus_prev", "unc_plus_recs")))


ggplot(tdf, aes(visit_numbers, diversity, col = method)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~uptake)

new_locs %>% 
  na.omit() %>% 
  mutate(locid = paste(lon, lat, sep = '_')) %>% 
  group_by(method, community, uptake, locid) %>% 
  summarise(sp = length(unique(species))) %>% 
  group_by(method, community, uptake, sp) %>% 
  tally()
