

library(tidyverse)
library(ggridges)
library(patchwork)

write = FALSE


## set percentage change
p_c <- 5


## sort data out
{
  ## load uptake values
  obvs_df <- read.csv('outputs/v4Community/v4n_new_obvs_perc_increase_1_50_spp50.csv', 
                      stringsAsFactors = FALSE)
  
  # load each of the evaluation files 
  cdf_0.1_uptake <- read.csv('outputs/v4Community/asv1_v4combined_outputs_comm1_50_spp50.csv', 
                             stringsAsFactors = FALSE) %>%
    mutate(uptake = '0.1',
           asv = 'asv1')
  cdf_0.01_uptake <- read.csv('outputs/v4Community/asv2_v4combined_outputs_comm1_50_spp50.csv', 
                              stringsAsFactors = FALSE) %>%
    mutate(uptake = '0.01',
           asv = 'asv2')
  
  cdf_0.5_uptake <- read.csv('outputs/v4Community/asv4_v4combined_outputs_comm1_50_spp50.csv', 
                             stringsAsFactors = FALSE) %>%
    mutate(uptake = '0.5',
           asv = 'asv4')
  
  cdf_0_uptake <- read.csv('outputs/v4Community/asv3_v4combined_outputs_comm1_50_spp50.csv', 
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
           delta_corr = corr - init_corr,
           delta_auc = auc - init_auc,
           delta_mean_sd = mean_sd - init_mean_sd,
           delta_max_sd = max_sd - init_max_sd)
  
  ## average across communities
  comm_df <- et %>%
    group_by(community, method, uptake) %>%
    summarise(mse = mean(mse, na.rm = T),
              corr = mean(corr, na.rm = T),
              auc = mean(auc, na.rm = T),
              mean_sd = mean(mean_sd, na.rm = T),
              max_sd = max(max_sd, na.rm = T),
              prev = median(prevalence, na.rm = T),
              init_mse = mean(init_mse, na.rm = T),
              init_corr = mean(init_corr, na.rm = T),
              init_auc = mean(init_auc, na.rm = T),
              init_mean_sd = mean(init_mean_sd, na.rm = T),
              init_max_sd = max(init_max_sd, na.rm = T),
              prev = median(prev)) %>%
    ungroup() %>% 
    mutate(delta_mse = sqrt(mse) - sqrt(init_mse),
           delta_corr = corr - init_corr,
           delta_auc = auc - init_auc,
           delta_mean_sd = mean_sd - init_mean_sd,
           delta_max_sd = max_sd - init_max_sd,
           method = factor(method, 
                           levels=c("initial", "none", "coverage", "prevalence",  
                                    "unc_plus_prev", "unc_plus_recs", "uncertainty")))
  
  
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
    summarise(n = sum(Observed, na.rm = T)) %>% 
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
           corr_cat = dplyr::ntile(init_corr, 10)) %>% 
    # ungroup() %>% 
    # rowwise() %>% 
    mutate(perc_inc_auc = (delta_auc)/(init_auc)*100,
           perc_inc_corr = (delta_corr)/(init_corr)*100,
           perc_inc_mse = (delta_mse)/(init_mse)*100,
           perc_imp_auc = ifelse(perc_inc_auc>= p_c, p_c, 
                                 ifelse(perc_inc_auc<= -p_c, -p_c, 0)),
           perc_imp_corr = ifelse(perc_inc_corr>=p_c, p_c, 
                                  ifelse(perc_inc_corr<= -p_c, -p_c, 0)),
           perc_imp_mse = ifelse(perc_inc_mse<= -p_c, p_c, 
                                 ifelse(perc_inc_mse>= p_c, -p_c, 0)),
           method = factor(method, 
                           levels=c("initial", "none", "coverage", "prevalence",  
                                    "unc_plus_prev", "unc_plus_recs", "uncertainty")))
  
  # number of models with > x% increase
  nmods <- etp %>%
    group_by(community, method, uptake) %>%
    summarise(n_mods_auc_1 = sum(perc_inc_auc>1, na.rm = T),
              n_mods_auc_2 = sum(perc_inc_auc>2, na.rm = T),
              n_mods_auc_5 = sum(perc_inc_auc>5, na.rm = T),
              n_mods_mse_1 = sum(perc_inc_mse< -1, na.rm = T),
              n_mods_mse_2 = sum(perc_inc_mse< -2, na.rm = T),
              n_mods_mse_5 = sum(perc_inc_mse< -5, na.rm = T),
              n_mods_corr_1 = sum(perc_inc_corr>1, na.rm = T),
              n_mods_corr_2 = sum(perc_inc_corr>2, na.rm = T),
              n_mods_corr_5 = sum(perc_inc_corr>5, na.rm = T)) 
  
  nmods_l <- pivot_longer(nmods, cols = 4:12) %>% 
    rowwise() %>% 
    mutate(eval_type = ifelse(grepl(x = name, pattern = 'auc'), 'auc',
                              ifelse(grepl(x = name, pattern = 'mse'), 'mse', 
                                     ifelse(grepl(x = name, pattern = 'corr'), 'corr', 'WRONG'))),
           inc_amount = as.numeric(gsub("[^\\d]+", "", name, perl=TRUE)),
           method = factor(method, 
                           levels=c("initial", "none", "coverage", "prevalence",  
                                    "unc_plus_prev", "unc_plus_recs", "uncertainty")))
  
}



## plotting -----


## figure 1 - model improvements + N models > %

cno1 <- ggplot(comm_df[comm_df$method!='initial' & comm_df$uptake!=0,] %>% 
                 mutate(facet = 'MSE'), 
               aes(x=method, y = delta_mse, fill = factor(uptake))) +
  geom_boxplot() +
  theme_bw() + 
  ylim(-0.008, 0.005) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  xlab('') +
  ylab('Delta MSE (lower = better)') +
  # theme(axis.text.x = element_blank()) +
  scale_fill_discrete(name = 'Uptake') +
  # facet_wrap(~facet) +
  scale_x_discrete(labels= c('Business\nas usual', 'Gap-filling', 'Rare species',
                             'Uncertainty of rare\nspecies', 'New areas\nand uncertainty',
                             'Uncertainty only')) +
  theme(text = element_text(size = 12))

cno1


msen <- ggplot(subset(nmods_l, uptake != 0 & inc_amount == 1 & method != 'initial' & eval_type == 'mse'), 
               aes(x = method, y = value, fill = factor(uptake))) +
  geom_boxplot() +
  scale_fill_discrete(name = "Uptake") +
  # facet_wrap(~eval_type,ncol = 3) +
  ylab('Number of models with\n>1% improvement') +
  xlab('') +
  # ggtitle('Number of models with 1% improvement for different uptake values') +
  theme_bw() +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  scale_x_discrete(labels= c('Business\nas usual', 'Gap-filling', 'Rare species',
                             'Uncertainty of rare\nspecies', 'New areas\nand uncertainty',
                             'Uncertainty only')) + 
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 25, vjust = 0.5))
msen


fig2 <- cno1 + msen + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'a')
fig2

if(write){
  ggsave(fig2, filename = 'outputs/plots/paper/figure2')
}




## figure 3 -- proportion of models

# change labels
meth_names <- list(
  "initial",
  "Business\nas usual",
  "Gap-filling",
  "Rare species",
  "Uncertainty of rare\nspecies",
  "New areas\nand uncertainty",
  "Uncertainty only"
)

etp_p <- subset(etp, method != 'initial')
levels(etp_p$method) <- unlist(meth_names)


ggplot(na.omit(subset(etp_p, uptake == 0.1 & method != 'initial')), 
       aes(x = prev_cat, fill = factor(perc_imp_mse))) +
  geom_bar(position="fill") +
  ylab('Proportion of models') +
  xlab('Prevalence category') +
  # ylim(0,0.25) +
  facet_grid(~method) +
  scale_fill_manual(name = 'Change in MSE', labels = c('< -5%', '0%', '> 5%'),
                    values = c("#E69F00", "#56B4E9", "#009E73")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1),
        text = element_text(size = 10))

# ggsave()
etp <- et %>% 
  # group_by(method, uptake) %>% 
  mutate(prev_cat = dplyr::ntile(prevalence, 10),
         auc_cat = dplyr::ntile(init_auc, 10),
         mse_cat = dplyr::ntile(init_mse, 10),
         corr_cat = dplyr::ntile(init_corr, 10)) %>% 
  # ungroup() %>% 
  # rowwise() %>% 
  mutate(perc_inc_auc = (delta_auc)/(init_auc)*100,
         perc_inc_corr = (delta_corr)/(init_corr)*100,
         perc_inc_mse = (delta_mse)/(init_mse)*100,
         perc_imp_auc = ifelse(perc_inc_auc>= p_c, p_c, 
                               ifelse(perc_inc_auc<= -p_c, -p_c, 0)),
         perc_imp_corr = ifelse(perc_inc_corr>=p_c, p_c, 
                                ifelse(perc_inc_corr<= -p_c, -p_c, 0)),
         perc_imp_mse = ifelse(perc_inc_mse<= -p_c, p_c, 
                               ifelse(perc_inc_mse>= p_c, -p_c, 
                                      ifelse(perc_inc_mse>= -p_c & perc_inc_mse< -1, -2.5,
                                             ifelse(perc_inc_mse<= p_c & perc_inc_mse> 1, 2.5, 0)))),
         method = factor(method, 
                         levels=c("initial", "none", "coverage", "prevalence",  
                                  "unc_plus_prev", "unc_plus_recs", "uncertainty")))

etp_p <- subset(etp, method != 'initial')
levels(etp_p$method) <- unlist(meth_names)


ggplot(na.omit(subset(etp_p, uptake == 0.1 & method != 'initial')), 
       aes(x = prev_cat, fill = factor(perc_imp_mse))) +
  geom_bar(position="fill") +
  ylab('Proportion of models') +
  xlab('Prevalence category') +
  # ylim(0,0.25) +
  facet_grid(~method) +
  scale_fill_manual(name = 'Change in MSE (%)', labels = c('< -5', '-5 to -1', 
                                                       '-1 to 1', '1 to 5', '> 5'),
                    values = c("#E69F00", "#56B4E9", "#009E73",
                               "#0072B2", "#D55E00")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1),
        text = element_text(size = 10))

