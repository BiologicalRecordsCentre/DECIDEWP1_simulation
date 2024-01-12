
## Script for the methodological figures 1 + 2

## packages
library(sf)
library(stars)
library(rnaturalearth)

library(tidyverse)
library(viridis)

# write figures to file?
write = TRUE

# methods used
method <- c("none", "uncertainty", "prevalence", "unc_plus_prev", "unc_plus_recs", "coverage")

#### data preparation
{
  # load the probability sampling layers
  {
    
    source('scripts/functions/filter_distance.R') ## crop within certain distance of centre
    
    # get GB map
    world <- ne_countries(scale = "medium", type = 'map_units', returnclass = "sf")
    uk <- world[world$admin == 'United Kingdom',]
    gb <- uk[uk$name != 'N. Ireland',] %>% 
      st_transform(crs = 27700)
    
    
    ### get cell weights (score for each cell)
    # list files
    fls <- list.files('outputs/comm4_asv1_investigating/', 
                      pattern = 'cellweights', full.names = TRUE)
    
    fls_names <- list.files('outputs/comm4_asv1_investigating/', 
                            pattern = 'cellweights')
    
    # load them all
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
    
    # convert to raster
    cwt_rast <- list()
    for(md in seq_along(method_unique)) {
      
      cell_wt_method <- cell_weights[cell_weights$method == method_unique[md],1:3]
      cwt_rast[[md]] <- raster::rasterFromXYZ(cell_wt_method)
      
    }
    
    cwt_rast[[1]] <- raster::crop(cwt_rast[[1]], cwt_rast[[2]])  
    
    all_comb <- raster::stack(cwt_rast)
    names(all_comb) <- method_unique  
    
    # crop
    crp_rst <- filter_distance(all_comb,
                               location = c(-2.2, 53.8),
                               distance = 20000)
    
    as_layers <- as.data.frame(crp_rst, xy = TRUE) %>% 
      pivot_longer(cols = 3:8)
    
  }
  
  # load the observation data
  {
    
    # get observations for each method
    method_locs <- data.frame()
    
    for(com in 1) { # community - only one needed
      
      for(a in 1){ # adaptive sampling version - only one needed
        
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
    
    # get community 1
    meth_locsC1 <- subset(method_locs, community == 1)
    
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
    
    # find the new locs that weren't the initial locs
    mlocsC1 <- meth_locsC1[!meth_locsC1$id %in% obvsinit$id,] %>% 
      dplyr::select(!species) %>% 
      distinct()
    
    # comine original and new locs
    mlocs_all <- rbind(mlocsC1, obvsinit)
    
    # get a subset
    mlocs_all_sub <- mlocs_all[!is.na(raster::extract(crp_rst[[2]], mlocs_all[,c(1,2)])), ]
    
  }
  
  # get subset of environmental data
  {
    env_data <- raster::stack('data/environmental_data/envdata_1km_no_corr_noNA.grd')[[22:24]]
    edf <- as.data.frame(env_data, xy = TRUE)
  }
  
  # get species and model data
  {
    # get community data
    c4 <- readRDS('outputs/v4Community/v4community_1_50_sim_initial.rds')
    c4cov <- readRDS('outputs/comm4_asv1_investigating/asv1_v4community_1_50_sim_initial_AS_uncertainty.rds')
    
    ## true prob of occurrence 
    c4sp3 <- c4[[5]]$pres_abs %>% 
      as.data.frame(xy=TRUE)
    c4covsp3 <- c4cov[[5]]$pres_abs %>% 
      as.data.frame(xy=TRUE)
    
    # model outputs
    # species 5
    load('outputs/comm4_asv1_investigating/asv1_v4lr_SDMs_GBnew_Sp5_initial_AS_uncertainty.rdata')
    pred_sp5 <- model_output$predictions
    # species 50
    load('outputs/comm4_asv1_investigating/asv1_v4lr_SDMs_GBnew_Sp50_initial_AS_uncertainty.rdata')
    pred_sp50 <- model_output$predictions
  }
}


#### figure 1 - plots for the AS methodological figure
{
  ### observations
  # sample to mimic only a few species
  uncert <- all_comb[['uncertainty']] %>% 
    as.data.frame(xy = TRUE) %>% 
    mutate(uncert = ifelse(is.na(uncertainty), NA, 1))
  
  initobs <- obvsinit #subset(mlocs_all, asv == 'init')
  
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
    ggsave(p_obs, filename = paste0('outputs/plots/paper/methods_plots/figure1_simmethodplot_sp_obs.png'),
           bg = 'transparent', width = 4.46, height = 7)
  }
  
  
  ### environmental layers
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
  
  # True species distributions
  true_sp <- ggplot() +
    geom_sf(data = gb, fill = 'white', col = 'white', size = 5) +
    geom_raster(data = c4sp3, aes(x,y,fill = layer)) +
    scale_fill_manual(na.value = 'transparent', values = c('grey', 'green')) +
    theme_void() +
    theme(legend.position = 'none')
  
  # species 5 distribution
  pred_sp_pl5 <- ggplot() +
    geom_sf(data = gb, fill = 'white', col = 'white', size = 4) +
    geom_raster(data = pred_sp5, aes(x,y,fill = mean)) +
    scale_fill_viridis(na.value = 'transparent') +
    theme_void() +
    theme(legend.position = 'none')
  
  # species 50 distribution
  pred_sp_pl50 <- ggplot() +
    geom_sf(data = gb, fill = 'white', col = 'white', size = 4) +
    geom_raster(data = pred_sp50, aes(x,y,fill = mean)) +
    scale_fill_viridis(na.value = 'transparent') +
    theme_void() +
    theme(legend.position = 'none')
  
  
  if(write) {  
    # environmental layers
    ggsave(ed1, filename = paste0('outputs/plots/paper/methods_plots/figure1_simmethodplot_antemp.png'),
           bg = 'transparent', width = 4.46, height = 7)  
    ggsave(ed2, filename = paste0('outputs/plots/paper/methods_plots/figure1_simmethodplot_mnditemp.png'),
           bg = 'transparent', width = 4.46, height = 7)  
    ggsave(ed3, filename = paste0('outputs/plots/paper/methods_plots/figure1_simmethodplot_isotherm.png'),
           bg = 'transparent', width = 4.46, height = 7)
    
    # species layers
    ggsave(true_sp, filename = paste0('outputs/plots/paper/methods_plots/figure1_simmethodplot_truepres.png'),
           bg = 'transparent', width = 4.46, height = 7)  
    
    ggsave(pred_sp_pl5, filename = paste0('outputs/plots/paper/methods_plots/figure1_simmethodplot_predpres_sp5.png'),
           bg = 'transparent', width = 4.46, height = 7)  
    
    ggsave(pred_sp_pl50, filename = paste0('outputs/plots/paper/methods_plots/figure1_simmethodplot_predpres_sp50.png'),
           bg = 'transparent', width = 4.46, height = 7)  
  }
}



#### figure 2 - plots for the AS methodological figure
{
  as_method <- method[2:5]
  
  for(md in 1:length(as_method)) {
    
    print(as_method[md])
    
    p <- ggplot(subset(as_layers, name == as_method[md]), 
                aes(x,y,fill = log(value))) +
      geom_raster()+
      scale_fill_viridis(na.value="transparent") +
      theme_void() +
      theme(legend.position = 'none') 
    
    if(write){
      ggsave(p, filename = paste0('outputs/plots/paper/methods_plots/figure2_ASmethodplot_', as_method[md], '.png'),
             bg = 'transparent', width = 2.2, height = 2.2)
    }
  }
  
  for(m in method[c(1,6)]){
    p_nas <- ggplot() +
      geom_raster(data = subset(as_layers, name == 'coverage'), aes(x,y,fill = log(value))) +
      geom_point(data = subset(mlocs_all_sub, method == 'initial'), aes(lon,lat), colour = 'black') +
      geom_point(data = subset(mlocs_all_sub, method == m), aes(lon,lat), colour = 'red') +
      scale_fill_viridis(na.value = 'transparent') +
      theme_void() +
      theme(legend.position = 'none')
    
    if(write){
      ggsave(p_nas, filename = paste0('outputs/plots/paper/methods_plots/figure2_ASmethodplot_', m, '.png'),
             bg = 'transparent', width = 2.2, height = 2.2)
    }
  }
}

