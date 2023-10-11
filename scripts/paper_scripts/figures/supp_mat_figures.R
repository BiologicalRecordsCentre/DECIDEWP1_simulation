
## Supplementary material figures

library(tidyverse)

write = TRUE

# Full method names for renaming
meth_names <- list(
  "initial",
  "Business\nas usual",
  "Gap-filling",
  "Rare species",
  "Uncertainty only",
  "Uncertainty of\nrare species",
  "Gap-filling\nwith uncertainty"
)

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
    ggsave(rankabund_butterfly, filename = 'outputs/plots/paper/supplementary_figures/supp_matt_rankabund_butterfly.png',
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
    arrange(X) %>%
    ungroup()
  
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
    ggsave(rank_abund_comm, filename = 'outputs/plots/paper/supplementary_figures/supp_matt_rankabund_communities_all.png',
           width = 9, height = 5)
  }  
}

#### Extra evaluation metrics

# create evaluation data frames
{
  
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
  
  # bind initial values to full dataset
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
  
  
  ## calculate differences between modelled values and initial data
  et <- et %>% 
    mutate(delta_mse = mse - init_mse,
           delta_medse = medianse - init_medse,
           delta_corr = corr - init_corr,
           delta_auc = auc - init_auc,
           delta_mean_sd = mean_sd - init_mean_sd,
           delta_max_sd = max_sd - init_max_sd)
  
  ## calculate average across communities
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
              prev = median(prev),
              
              # mean of deltas
              delta_mse = mean(delta_mse, na.rm = TRUE),
              delta_medse = mean(delta_mse, na.rm = TRUE),
              delta_corr = mean(delta_corr, na.rm = TRUE),
              delta_auc = mean(delta_auc, na.rm = TRUE),
              delta_mean_sd = mean(delta_mean_sd, na.rm = TRUE),
              delta_max_sd = mean(delta_max_sd, na.rm = TRUE)) %>%
    ungroup() %>% 
    mutate(method = factor(method, 
                           levels=c("initial", "none", "coverage", "prevalence",  
                                    "uncertainty", "unc_plus_prev", "unc_plus_recs")))
  
  
  
  ## get model improvements by certain percentages
  # percentage 
  p_c = 5
  
  
  etp <- et %>% 
    # split data into different categories
    mutate(prev_cat = as.numeric(cut_number(prevalence,10)), #dplyr::ntile(prevalence, 10), # prevalence
           auc_cat = as.numeric(cut_number(init_auc,10)), #dplyr::ntile(init_auc, 10),
           mse_cat = as.numeric(cut_number(init_mse,10)), #dplyr::ntile(init_mse, 10),
           medse_cat = as.numeric(cut_number(init_medse,10)), #dplyr::ntile(init_medse, 10),
           corr_cat = as.numeric(cut_number(init_corr,10)), #dplyr::ntile(init_corr, 10))
           
           # get percentage increase
           perc_inc_auc = (delta_auc)/(init_auc)*100,
           perc_inc_corr = (delta_corr)/(init_corr)*100,
           perc_inc_mse = (delta_mse)/(init_mse)*100,
           perc_inc_medse = (delta_medse)/(init_medse)*100,
           
           # Get number of models per percentage increase increment
           perc_imp_auc = ifelse(perc_inc_auc<= -p_c, -p_c, 
                                 ifelse(perc_inc_auc>= p_c, p_c, 
                                        ifelse(perc_inc_auc>= -p_c & perc_inc_auc< 0, -2.5,
                                               ifelse(perc_inc_auc<= p_c & perc_inc_auc> 0, 2.5, 0)))),
           perc_imp_corr = ifelse(perc_inc_corr<= -p_c, -p_c, 
                                  ifelse(perc_inc_corr>= p_c, p_c, 
                                         ifelse(perc_inc_corr>= -p_c & perc_inc_corr< 0, -2.5,
                                                ifelse(perc_inc_corr<= p_c & perc_inc_corr> 0, 2.5, 0)))),
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
                                    "uncertainty", "unc_plus_prev", "unc_plus_recs"))) %>%
    ungroup()
  
  # remove initial data
  etp_p2 <- subset(etp, method != 'initial')
  levels(etp_p2$method) <- unlist(meth_names)
  
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
              n_mods_corr_5 = sum(perc_inc_corr>5, na.rm = TRUE)) %>%
    ungroup()
  
  # pivot to long format
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
  
}

# create proportion increase of range data frame
{
  
  ## load all the observations for each of the AS versions
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
  
  # bind all together, reorder and rename the methods
  df <- rbind(asv1,asv2,asv3,asv4) %>% 
    mutate(method = factor(method, 
                           levels=c("initial", "none", "coverage", "prevalence",  
                                    "uncertainty", "unc_plus_prev", "unc_plus_recs")))
  levels(df$method) <- unlist(meth_names)
  
  
  # extract the initial locations
  init_nona <- df[which(df$method == 'initial'),]
  
  # remove all initial locations so that we're only looking at new locations
  new_locs <- df[!df$id %in% init_nona$id,]
  
  
  # Calculate Proportion of a species' total range covered by new locations from each sampling method
  ## How much of the species' true range is sampled using each of the sampling methods?
  # better to look at the change because that will be dwarfed by the total number of observations
  # this is just the number of new observations that I've already calculated in new_locs
  
  
  ### first, get the number of cells that every species occurs in
  # use the prevalence * n_cells_uk
  
  ## to get the total number of cells in the uk
  # get environmental data
  library(terra)
  
  # climate layers cover all cells
  env_data <- rast('data/environmental_data/envdata_1km_no_corr_noNA.grd')[[22]] 
  edf <- as.data.frame(env_data, xy=TRUE)
  head(edf)
  dim(edf) # nrows = n_cells that aren't NA
  
  # calculate number of cells that each species occurs in
  new_locs$prevalence_ncells <- new_locs$prevalence * dim(edf)[1]
  
  # calculate the proportion range that new observations cover
  prop_range <- new_locs %>%
    mutate(prev_cat = as.numeric(cut_number(prevalence,10))) %>% # prevalence category)
    na.omit() %>% 
    group_by(community, uptake, method, species, prevalence_ncells) %>% 
    summarise(new_obs = sum(Observed), # total number of new observations
              prevalence = unique(prevalence), # original prevalence
              prev_cat = unique(prev_cat),
              prevalence_ncells = unique(prevalence_ncells))  %>%
    ungroup() %>% 
    mutate(prop_cov_increase = new_obs/prevalence_ncells, # change in the proportion of the TRUE range that is sampled
           prop_cov_inc_cat = dplyr::ntile(prop_cov_increase, 4),
           sp_id = paste(community, uptake, species, sep = "_"))
  
}


#### Plotting -----

## figure 3 - model improvements + N models > %
{
  # MSE with outliers
  {
    ## go through all model evaluation methods
    # figure 3a - model improvements
    msewithoutliers <- ggplot(comm_df[comm_df$method!='initial' & comm_df$uptake!=0,] %>% 
                                mutate(facet = 'MSE'), 
                              aes(x=method, y = delta_mse, fill = factor(uptake))) +
      geom_boxplot() +
      theme_bw() + 
      # ylim(-(2*sd(comm_df$delta_mse)), 2*sd(comm_df$delta_mse)) +
      geom_hline(yintercept = 0, linetype = 'dashed') +
      xlab('') +
      ylab('Delta MSE (lower = better)') +
      scale_fill_manual(name = 'Uptake (%)', labels = c(1, 10, 50),
                        values = c("#E69F00", "#56B4E9", "#009E73")) +
      scale_x_discrete(labels= c('Business \n as usual', 'Gap-filling', 'Rare species',
                                 'Uncertainty only', 'Uncertainty of \n rare species', 
                                 'Gap-filling \n with uncertainty')) +
      theme(text = element_text(size = 12))
    msewithoutliers
    
  }
  
  # Correlation
  {
    # figure 3a - model improvements
    corrimp <- ggplot(comm_df[comm_df$method!='initial' & comm_df$uptake!=0,] %>% 
                        mutate(facet = 'Correlation'), 
                      aes(x=method, y = delta_corr, fill = factor(uptake))) +
      geom_boxplot() +
      theme_bw() + 
      geom_hline(yintercept = 0, linetype = 'dashed') +
      xlab('') +
      ylab('Delta Correlation (higher = better)') +
      scale_fill_manual(name = 'Uptake (%)', labels = c(1, 10, 50),
                        values = c("#E69F00", "#56B4E9", "#009E73")) +
      scale_x_discrete(labels= c('Business \n as usual', 'Gap-filling', 'Rare species',
                                 'Uncertainty only', 'Uncertainty of \n rare species', 
                                 'Gap-filling \n with uncertainty')) +
      theme(text = element_text(size = 12),
            axis.text.x = element_blank())
    
    corrimp
    
    # figure 3b - number of models >X%
    corrperc <- ggplot(subset(nmods_l, uptake != 0 & inc_amount == 1 & method != 'initial' & eval_type == 'corr'), 
                       aes(x = method, y = value, fill = factor(uptake))) +
      geom_boxplot() +
      scale_fill_manual(name = "Uptake (%)", labels = c(1, 10, 50),
                        values = c("#E69F00", "#56B4E9", "#009E73")) +
      ylab('Number of models with\n>1% improvement') +
      xlab('') +
      theme_bw() +
      scale_x_discrete(labels= c('Business \n as usual', 'Gap-filling', 'Rare species',
                                 'Uncertainty\nonly', 'Uncertainty of \n rare species', 
                                 'Gap-filling \n with uncertainty')) + 
      theme(text = element_text(size = 12),
            axis.text.x = element_text(size = 12, angle = 0, vjust = 0))
    corrperc
    
    # combine them all
    fig3corr <- corrimp / corrperc + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'a')
    fig3corr
    
  }
  
  # AUC
  {
    # figure 3a - model improvements
    aucimp <- ggplot(comm_df[comm_df$method!='initial' & comm_df$uptake!=0,] %>% 
                       mutate(facet = 'AUC'), 
                     aes(x=method, y = delta_auc, fill = factor(uptake))) +
      geom_boxplot() +
      theme_bw() + 
      geom_hline(yintercept = 0, linetype = 'dashed') +
      xlab('') +
      ylab('Delta AUC (higher = better)') +
      scale_fill_manual(name = 'Uptake (%)', labels = c(1, 10, 50),
                        values = c("#E69F00", "#56B4E9", "#009E73")) +
      scale_x_discrete(labels= c('Business \n as usual', 'Gap-filling', 'Rare species',
                                 'Uncertainty only', 'Uncertainty of \n rare species', 
                                 'Gap-filling \n with uncertainty')) +
      theme(text = element_text(size = 12),
            axis.text.x = element_blank())
    
    aucimp
    
    # figure 3b - number of models >X%
    aucperc <- ggplot(subset(nmods_l, uptake != 0 & inc_amount == 1 & method != 'initial' & eval_type == 'auc'), 
                      aes(x = method, y = value, fill = factor(uptake))) +
      geom_boxplot() +
      scale_fill_manual(name = "Uptake (%)", labels = c(1, 10, 50),
                        values = c("#E69F00", "#56B4E9", "#009E73")) +
      ylab('Number of models with\n>1% improvement') +
      xlab('') +
      theme_bw() +
      scale_x_discrete(labels= c('Business \n as usual', 'Gap-filling', 'Rare species',
                                 'Uncertainty\nonly', 'Uncertainty of \n rare species', 
                                 'Gap-filling \n with uncertainty')) + 
      theme(text = element_text(size = 12),
            axis.text.x = element_text(size = 12, angle = 0, vjust = 0))
    aucperc
    
    # combine them all
    fig3auc <- aucimp / aucperc + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'a')
    fig3auc
    
  }
  
  if(write){
    ggsave(msewithoutliers, filename = 'outputs/plots/paper/supplementary_figures/figure3_mse_outliers.png',
           width = 9.5, height = 4.5)
    ggsave(fig3corr, filename = 'outputs/plots/paper/supplementary_figures/figure3_improve_uptake_corr.png',
           width = 9.5, height = 9.5)
    ggsave(fig3auc, filename = 'outputs/plots/paper/supplementary_figures/figure3_improve_uptake_auc.png',
           width = 9.5, height = 9.5)
  }
  
}


## figure 4 - proportion of models in different percentage categories
{
  # mse
  {
    fig4_mse_uptake <- ggplot(na.omit(subset(etp_p2, uptake != 0 & method != 'initial')), 
                              aes(x = prev_cat, fill = factor(perc_imp_mse))) +
      geom_bar(position="fill") +
      ylab('Proportion of models') +
      xlab('Prevalence category') +
      facet_grid(uptake~method) +
      scale_fill_manual(name = 'Improvement in\nmodel MSE (%)', 
                        labels = c('< -5', '-5 to 0', 
                                   '0 to 5', '> 5'),
                        values = c("#E69F00", "#56B4E9", "#009E73",
                                   "#0072B2")) +#, "#D55E00")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1),
            text = element_text(size = 14))
    fig4_mse_uptake
  }
  
  # correlation
  {
    fig4_corr_uptake <- ggplot(na.omit(subset(etp_p2, uptake != 0 & method != 'initial')), 
                               aes(x = prev_cat, fill = factor(perc_imp_corr))) +
      geom_bar(position="fill") +
      ylab('Proportion of models') +
      xlab('Prevalence category') +
      facet_grid(uptake~method) +
      scale_fill_manual(name = 'Improvement in\nmodel correlation (%)', 
                        labels = c('< -5', '-5 to 0', 
                                   '0 to 5', '> 5'),
                        values = c("#E69F00", "#56B4E9", "#009E73",
                                   "#0072B2")) +#, "#D55E00")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1),
            text = element_text(size = 14))
    fig4_corr_uptake
  }
  
  # auc
  {
    fig4_auc_uptake <- ggplot(na.omit(subset(etp_p2, uptake != 0 & method != 'initial')), 
                              aes(x = prev_cat, fill = factor(perc_imp_auc))) +
      geom_bar(position="fill") +
      ylab('Proportion of models') +
      xlab('Prevalence category') +
      facet_grid(uptake~method) +
      scale_fill_manual(name = 'Improvement in\nmodel AUC (%)', 
                        labels = c('< -5', '-5 to 0', 
                                   '0 to 5', '> 5'),
                        values = c("#E69F00", "#56B4E9", "#009E73",
                                   "#0072B2")) +#, "#D55E00")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1),
            text = element_text(size = 14))
    fig4_auc_uptake
  }
  
  if(write){
    ggsave(fig4_mse_uptake, filename = 'outputs/plots/paper/supplementary_figures/fig4_mse_uptake.png',
           width = 11, height = 10)
    ggsave(fig4_corr_uptake, filename = 'outputs/plots/paper/supplementary_figures/fig4_corr_uptake.png',
           width = 11, height = 10)
    ggsave(fig4_auc_uptake, filename = 'outputs/plots/paper/supplementary_figures/fig4_auc_uptake.png',
           width = 11, height = 10)
  }
}


## figure 5 - Proportion of a species' total range covered by new locations from each sampling method
{
  ## How much of the species' true range is sampled using each of the sampling methods?
  fig5_proprange_full <- ggplot(subset(prop_range, uptake!=0), aes(x = method, y = prop_cov_increase, fill = factor(uptake))) +
    geom_boxplot() +
    # ylim(0, 0.04)+ #3*sd(prop_range$prop_cov_increase)) +
    ylab("Increase in coverage of species ranges\nfrom new observations") +
    xlab('') +
    scale_fill_manual(name = 'Uptake (%)', labels = c(1, 10, 50),
                      values = c("#E69F00", "#56B4E9", "#009E73")) +
    theme_classic() +
    theme(text = element_text(size = 12))
  fig5_proprange_full
  
  ## How much of the species' true range is sampled using each of the sampling methods?
  fig5_proprange_3sd <- ggplot(subset(prop_range, uptake!=0), aes(x = method, y = prop_cov_increase, fill = factor(uptake))) +
    geom_boxplot() +
    ylim(0, 3*sd(prop_range$prop_cov_increase)) +
    ylab("Increase in coverage of species ranges\nfrom new observations") +
    xlab('') +
    scale_fill_manual(name = 'Uptake (%)', labels = c(1, 10, 50),
                      values = c("#E69F00", "#56B4E9", "#009E73")) +
    theme_classic() +
    theme(text = element_text(size = 12))
  fig5_proprange_3sd
  
  fig5_proprange <- fig5_proprange_full/fig5_proprange_3sd + 
    plot_annotation(tag_levels = 'a') + 
    plot_layout(guides = 'collect')
  fig5_proprange
  
  
  # Proportion of range increase by prevalence
  prop_range %>% 
    group_by(prop_cov_inc_cat) %>% 
    reframe(pinc = range(prop_cov_increase)) %>% 
    mutate(pinc*100) #### Create custom groups using this?
  
  prop_range_prevpl <- ggplot(na.omit(subset(prop_range, uptake != 0 & method != 'initial')), 
                              aes(x = prev_cat, fill = factor(prop_cov_inc_cat))) +
    geom_bar(position="fill") +
    ylab('Proportion of species') +
    xlab('Prevalence category') +
    facet_grid(uptake~method) +
    scale_fill_manual(name = 'Proportion of range sampled\nby new observations (%)',
                      labels = c('0.001 < prop > 0.01', 
                                 '0.01 < prop > 0.03',
                                 '0.03 < prop > 0.05',
                                 '0.05 < prop > 6'),
                      values = c("#E69F00", "#56B4E9", "#009E73",
                                 "#0072B2")) +#, "#D55E00")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1),
          text = element_text(size = 14))
  prop_range_prevpl
  
  if(write){
    ggsave(fig5_proprange, filename = 'outputs/plots/paper/supplementary_figures/fig5_proprange.png',
           width = 9.5, height = 9.5)
    ggsave(prop_range_prevpl, filename = 'outputs/plots/paper/supplementary_figures/fig5_proprange_prevalence.png',
           width = 11.5, height = 8.5)
  }
}

## Raw AUC/MSE/correlations scores for models before/after AS? 
head(et)
eval_comp <- subset(et, method != 'initial') %>% 
  dplyr::select(method, community, uptake, species, mse, corr, auc, init_mse, init_auc, init_corr)
head(eval_comp)

auc <- eval_comp %>% 
  dplyr::select(method, community, uptake, auc, init_auc) %>% 
  pivot_longer(cols = c(auc, init_auc))

head(auc)

ggplot(auc, aes(method, value, fill = name)) +
  geom_boxplot() +
  facet_wrap(~uptake)

