library(tidyverse)

c <- readRDS('outputs/v4Community/asv1_v4community_1_50_sim_initial_AS_coverage.rds')
c[[1]]$model_variables


cn <- readRDS('outputs/v4Community/asv1_v4community_1_50_sim_initial_AS_none.rds')
cn[[1]]$model_variables

ci <- readRDS('outputs/v4Community/v4community_1_50_sim_initial.rds')
ci[[1]]$model_variables



# load each of the evaluation files 
cdf_0.1_uptake <- read.csv('outputs/v4Community/asv1_v4combined_outputs_comm1_50_spp50.csv', 
                           stringsAsFactors = FALSE) %>%
  mutate(uptake = '0.1')
cdf_0.01_uptake <- read.csv('outputs/v4Community/asv2_v4combined_outputs_comm1_50_spp50.csv', 
                            stringsAsFactors = FALSE) %>%
  mutate(uptake = '0.01')

cdf_0.5_uptake <- read.csv('outputs/v4Community/asv4_v4combined_outputs_comm1_50_spp50.csv', 
                           stringsAsFactors = FALSE) %>%
  mutate(uptake = '0.5')

cdf_0_uptake <- read.csv('outputs/v4Community/asv3_v4combined_outputs_comm1_50_spp50.csv', 
                         stringsAsFactors = FALSE) %>%
  mutate(uptake = '0')

# combine all the files
cdf <- rbind(cdf_0_uptake, cdf_0.1_uptake,cdf_0.01_uptake,cdf_0.5_uptake)



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
up <- et %>% 
  mutate(delta_mse = mse - init_mse,
         delta_corr = corr - init_corr,
         delta_auc = auc - init_auc,
         delta_mean_sd = mean_sd - init_mean_sd,
         delta_max_sd = max_sd - init_max_sd) #%>% 
# filter(uptake == '0.1') 
head(up)


outt <- data.frame()
outdf <- data.frame()

for(i in 1:50){
  t <- up[up$species == paste0('Sp', i),]
  # head(t)
  
  t$vars <- paste(ci[[i]]$model_variables, collapse = ', ')
  # head(t)
  
  t2 <- t %>%
    mutate(col1 = strsplit(as.character(vars), ", ")) %>%
    unnest(col1) %>%
    filter(col1 != "")
  
  outt <- rbind(outt, t)
  
  outdf <- rbind(outdf, t2)
}

outt <- outt %>% 
  group_by(uptake, method) %>% 
  mutate(med = median(delta_mse, na.rm = T)) %>% 
  ungroup()

ggplot(outt, aes(x = reorder(vars, med), y = delta_mse)) +
  geom_boxplot()

head(outdf)

cats <- data.frame(vars = unique(outdf$col1), 
                   type = c('hab','hab','hab',
                            'hab','clim','clim',
                            'hab','hab','hab',
                            'hab','hab','hab',
                            'terr','clim','hab',
                            'hab','clim','clim',
                            'clim','hab','terr',
                            'hab','clim','hab',
                            'hab','terr','hab',
                            'hab','hab','hab',
                            'clim','clim','hab'))
head(cats)

# outdf <- 
outdf <- outdf %>% 
  rowwise() %>% 
  mutate(type = cats$type[col1 == cats$vars])


ggplot(outdf[outdf$method != 'initial',], aes(x = method, y = delta_mse, fill = col1)) +
  geom_boxplot() + #outlier.shape = NA) +
  # ylim(-0.01, 0.0175) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  facet_wrap(~uptake)

ggplot(outdf[outdf$method != 'initial',], aes(x = method, y = delta_mse, fill = type)) +
  geom_boxplot() + #outlier.shape = NA) +
  # ylim(-0.01, 0.0175) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  facet_wrap(~uptake)

ggplot(outdf[outdf$method != 'initial' & outdf$uptake == '0.01',], aes(x = method, y = delta_mse, fill = type)) +
  geom_boxplot() + #outlier.shape = NA) +
  # ylim(-0.05, 0.05) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  facet_wrap(~uptake)





