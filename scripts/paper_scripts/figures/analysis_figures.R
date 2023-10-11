
## Script for analysis figures 3,4 and 5


library(tidyverse)
library(patchwork)

#### Data preparation

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

write = TRUE

## create evaluation data frames
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
                                    "uncertainty", "unc_plus_prev", "unc_plus_recs")))
  
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
              n_mods_corr_5 = sum(perc_inc_corr>5, na.rm = TRUE)) 
  
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


## create observations data frame
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
  
  # calculate number of observations in new locations for each community
  comm_new_locs <- new_locs %>% 
    mutate(locid = paste(lon, lat, sep = '_')) %>% 
    group_by(method, community, uptake) %>% 
    summarise(unique_new_locs_community = length(unique(locid)))
  
  # calculate the number of observations in new locations for each species
  sp_new_locs <- new_locs %>% 
    mutate(locid = paste(lon, lat, sep = '_')) %>% 
    group_by(method, community, species, prevalence, asv, uptake) %>% 
    summarise(new_locations_species = length(unique(locid)))
  
}


#### figure 3 - model improvements + N models > %
{
  
  ## go through all model evaluation methods
  # figure 3a - model improvements
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
    scale_x_discrete(labels= c('Business \n as usual', 'Gap-filling', 'Rare species',
                               'Uncertainty only', 'Uncertainty of \n rare species', 
                               'Gap-filling \n with uncertainty')) +
    theme(text = element_text(size = 12),
          axis.text.x = element_blank())
  
  cno1
  
  # figure 3b - number of models >X%
  msen <- ggplot(subset(nmods_l, uptake != 0 & inc_amount == 1 & method != 'initial' & eval_type == 'mse'), 
                 aes(x = method, y = value, fill = factor(uptake))) +
    geom_boxplot() +
    scale_fill_manual(name = "Uptake (%)", labels = c(1, 10, 50),
                      values = c("#E69F00", "#56B4E9", "#009E73")) +
    # facet_wrap(~eval_type,ncol = 3) +
    ylab('Number of models with\n>1% improvement') +
    xlab('') +
    theme_bw() +
    scale_x_discrete(labels= c('Business \n as usual', 'Gap-filling', 'Rare species',
                               'Uncertainty\nonly', 'Uncertainty of \n rare species', 
                               'Gap-filling \n with uncertainty')) + 
    theme(text = element_text(size = 12),
          axis.text.x = element_text(size = 12, angle = 0, vjust = 0))
  msen
  
  # combine them all
  fig3 <- cno1 / msen + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'a')
  fig3
  
  if(write){
    ggsave(fig3, filename = 'outputs/plots/paper/analysis_plots/figure3_improve_uptake.png',
           width = 9.5, height = 9.5)
  }
}



#### figure 4 - proportion of models in different percentage categories
{
  fig4 <- ggplot(na.omit(subset(etp_p2, uptake == 0.5 & method != 'initial')), 
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
  fig4
  
  
  if(write){
    ggsave(fig4, filename = 'outputs/plots/paper/analysis_plots/figure4_propmods_0_1_5.png',
           width = 11, height = 5)
  }
  
}



#### figure 5 - number of new locations by community and by species
{
  # number of new locations for each species
  ## How many new locations are  each species seen in? (one observation in each location only though)
  fig5a <- ggplot(data = subset(sp_new_locs, method != 'initial' & uptake != 0), 
                  aes(x=method, y=new_locations_species, fill = factor(uptake))) +
    geom_boxplot() +
    ylab('Number of observations in\nnew locations per species') +
    xlab('') +
    scale_fill_manual(name = 'Uptake (%)', labels = c(1, 10, 50),
                      values = c("#E69F00", "#56B4E9", "#009E73")) +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          text = element_text(size = 12))
  
  # number of new locations per community
  ## in a community, how many new locations are visited?
  fig5b <- ggplot(data = subset(comm_new_locs, method != 'initial' & uptake != 0), 
                  aes(x=method, y=unique_new_locs_community, fill = factor(uptake))) +
    geom_boxplot() +
    ylab('Number of observations in\nnew locations per community') +
    xlab('') +
    scale_fill_manual(name = 'Uptake (%)', labels = c(1, 10, 50),
                      values = c("#E69F00", "#56B4E9", "#009E73")) +
    theme_classic() +
    theme(text = element_text(size = 12))
  
  
  fig5 <- fig5a/fig5b + 
    plot_annotation(tag_levels = 'a') + 
    plot_layout(guides = 'collect')
  fig5
  
  if(write){
    ggsave(fig5, filename = 'outputs/plots/paper/analysis_plots/figure5_new_locs.png',
           width = 9.5, height = 9.5)
  }
}








##




############# Ignore rest of this I think

# number of new sightings vs prevalence of the species
ggplot(data = subset(sp_new_locs, method != 'initial' & uptake == 0.5), 
       aes(x=prevalence, y=new_observations_species, colour = factor(uptake))) +
  geom_point() +
  facet_wrap(~factor(method)) +
  theme_bw()

# the prevalence of species at new sites
species_summs_df <- new_locs %>% 
  na.omit() %>% 
  mutate(locid = paste(lon, lat, sep = '_')) %>% 
  group_by(method, community, uptake, locid) %>% # group by site in each new 
  summarise(nsp = length(species),
            prev_per_loc = median(prevalence),
            unique_species_per_loc = length(unique(species)),
            num_species_per_loc = sum(Observed)) # %>% 
# mutate(method = factor(method, 
#                        levels=c("initial", "none", "coverage", "prevalence",  
#                                 "uncertainty", "unc_plus_prev", "unc_plus_recs")))

# species_summs_df <- species_summs_df[!duplicated(species_summs_df),]

ggplot(subset(species_summs_df, method != 'initial' & uptake != 0),
       aes(x = method, y = prev_per_loc, fill = factor(uptake))) +
  geom_boxplot()


ggplot(subset(species_summs_df, method != 'initial' & uptake != 0),
       aes(x = method, y = unique_species_per_loc, fill = factor(uptake))) +
  geom_boxplot()


ggplot(subset(species_summs_df, method != 'initial' & uptake != 0),
       aes(x = method, y = num_species_per_loc, fill = factor(uptake))) +
  geom_boxplot()




### community level stuff
loc_summ <- new_locs %>% 
  na.omit() %>% 
  mutate(locid = paste(lon, lat, sep = '_')) %>% 
  group_by(method, community, uptake) %>% # work out some community-level summaries
  summarise(unique_locs = length(unique(locid)), # number of new sites visited for each community
            total_obs = sum(Observed),
            diversity = length(unique(species))) #%>%  # n unique species
#   mutate(method = factor(method, 
#                          levels=c("initial", "none", "coverage", "prevalence",  
#                                   "uncertainty", "unc_plus_prev", "unc_plus_recs")))
# levels(loc_summ$method) <- unlist(meth_names)



# number of new sites visited in each community
ggplot(data = subset(loc_summ, method != 'initial' & uptake != 0), 
       aes(x=method, y=unique_locs, fill = factor(uptake))) +
  geom_boxplot() +
  ylab('Number of new locations visited\nper community') +
  xlab('') +
  scale_fill_manual(name = 'Uptake (%)', labels = c(1, 10, 50),
                    values = c("#E69F00", "#56B4E9", "#009E73")) +
  theme_classic() +
  theme(#axis.text.x = element_blank(),
    text = element_text(size = 12))



# total number of observations from new locations
ggplot(data = subset(loc_summ, method != 'initial' & uptake != 0), 
       aes(x=method, y=total_obs, fill = factor(uptake))) +
  geom_boxplot() +
  ylab('Total number of observations\nfrom new locations') +
  xlab('') +
  scale_fill_manual(name = 'Uptake (%)', labels = c(1, 10, 50),
                    values = c("#E69F00", "#56B4E9", "#009E73")) +
  theme_classic() +
  theme(#axis.text.x = element_blank(),
    text = element_text(size = 12))






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
  summarise(prev_per_loc = median(prevalence),
            unique_species_per_loc = length(unique(species)),
            num_species_per_loc = sum(Observed)) %>% 
  mutate(method = factor(method, 
                         levels=c("initial", "none", "coverage", "prevalence",  
                                  "uncertainty", "unc_plus_prev", "unc_plus_recs")))






ggplot(data = subset(sp_new_locs,method != 'initial' & uptake != 0), 
       aes(x=method, y=n, fill = factor(uptake))) +
  geom_boxplot() +
  ylab('Number of new locations visited') +
  xlab('') +
  theme_classic() +
  theme(#axis.text.x = element_blank(),
    text = element_text(size = 12))

## plots
# total number of observations per community
s1 <- ggplot(data = subset(loc_summ,method != 'initial' & uptake == 0.5), 
             aes(x=method, y=unique_locs)) +
  geom_boxplot() +
  ylab('Number of new locations visited') +
  xlab('') +
  theme_classic() +
  theme(#axis.text.x = element_blank(),
    text = element_text(size = 12))

s2 <- ggplot(data = subset(loc_summ, uptake==0.5), 
             aes(x=method, y=av_div_per_loc)) +
  geom_boxplot() +
  ylab('Species diversity\nper new location') +
  xlab('') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        text = element_text(size = 12))

s3 <- ggplot(data = subset(prev_df_unique, uptake==0.5), 
             aes(x=method, y=med_prev)) +
  geom_boxplot() +
  ylab('Median prevalence\nat new locations') +
  xlab('') +
  theme_classic() +
  theme(text = element_text(size = 12))

fig5 <- s1/s2/s3 +
  plot_annotation(tag_levels = 'a')
fig5
