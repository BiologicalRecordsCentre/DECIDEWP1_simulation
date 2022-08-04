## working out whether rare species are more likely to be seen with 
## certain AS methods
library(tidyverse)

## to be run on jasmin
for(i in 2:4){
  method = c("none", "uncertainty", "prevalence", "unc_plus_prev", "unc_plus_recs", "coverage")
  
  AS_version = paste0('asv',i)
  comm_version = 'v4'
  
  # get initial observations for each species
  
  comm_tally <- data.frame()
  
  for(c in 1:50){
    
    print(paste('community', c))
    
    fls_ls <- grep(AS_version, list.files(paste0('Outputs/', comm_version, 'communities_1km/', comm_version, 'community_', c, '_50_sim/preds_and_obsvs/'),
                                          pattern = '_sim_observations.csv', full.names = T), value = T)
    
    obvsinit <- readRDS(paste0('Outputs/v4communities_1km/', comm_version, 'community_', c, '_50_sim/', comm_version, 'community_',c,'_50_sim_initial.rds'))
    
    spp_obs <- data.frame()
    
    for(sp in 1:50){
      
      if(file.size(fls_ls[sp])==3) next() ## skip species when no observations - not ideal... 
      
      # initial observations
      init_obs <- obvsinit[[sp]]$observations %>% #[!is.na(obvsinit[[sp]]$observations$Observed),] %>% 
        mutate(method = 'initial',
               id = paste(lon, lat))
      # print(dim(init_obs))
      
      
      # all observations after AS
      o <- read.csv(fls_ls[sp]) %>% 
        mutate(id = paste(lon, lat)) %>% 
        filter(!is.na(Observed))
      # print(dim(o)) ## this is a bigger DF than initial obs because it has all methods  
      
      # remove initial observations
      new_locs <- o[!o$id %in% init_obs$id,] %>% 
        group_by(community, method, species, prevalence) %>%
        summarise(n_obs = sum(Observed, na.rm = T))
      
      spp_obs <- rbind(spp_obs, o)
      
      # print(new_locs)
    }
    
    comm_tally <- rbind(comm_tally, spp_obs)
    
    
  }
  
  write.csv(comm_tally, file = paste0('Outputs/', comm_version, 'communities_1km/', AS_version, '_', comm_version, 'all_observations.csv'))
  
}



obs <- read.csv('outputs/v4Community/asv1_v4n_observations_by_group.csv')
head(obs)


ggplot(obs, aes(x = prevalence, y = n_obs, colour = method)) +
  # geom_point()
  geom_smooth()

rare <- obs %>% 
  mutate(cat = cut(prevalence, breaks = quantile(prevalence)),
         prob = n_obs/2000)

ggplot(rare, aes(method, prob)) +
  geom_boxplot() +
  facet_wrap(~cat, scales = 'free')

ggplot(rare, aes(method, prob)) +
  geom_boxplot()

obs %>% 
  mutate(ratio = n_obs/prevalence) %>% 
  ggplot(aes(x=method, y = ratio)) +
  geom_boxplot()

qs <- obs %>% 
  group_by(method) %>% 
  summarise(lwr = quantile(n_obs)[2],
            upr = quantile(n_obs)[4])

ggplot(obs, aes(n_obs)) + 
  geom_histogram() +
  geom_vline(data = qs, aes(xintercept = lwr)) +
  geom_vline(data = qs, aes(xintercept = upr)) +
  facet_wrap(~method)


unique(obs$community)


##########Need to do a github issue about this!!!

### looking at the number of observations per method and prevalence
library(tidyverse)
library(patchwork)

df <- read_csv('outputs/v4Community/asv1_v4all_observations.csv')
head(df)

# initial locations
init <- df[df$method == 'initial',]

# all others
meths <- df[df$method != 'initial',]

# remove initial locations - only looking at new locations
new_locs <- meths[!meths$id %in% init$id,]

## some summaries
loc_summ <- new_locs %>% 
  na.omit() %>% 
  group_by(method, community, species, prevalence) %>% 
  summarise(n = sum(Observed)) %>% # number of observations of each species for all methods in each community
  group_by(method, community) %>% # work out some community-level summaries
  summarise(total_obs = sum(n),
            diversity = length(unique(species)), # n unique species
            prev = median(prevalence)) # median prevalence


## plots
s1 <- ggplot(data = loc_summ, aes(x=method, y=total_obs)) +
  geom_boxplot() +
  ylab('Total observations') +
  theme_classic()
s2 <- ggplot(data = loc_summ, aes(x=method, y=diversity)) +
  geom_boxplot() +
  ylab('Unique species') +
  theme_classic()
s3 <- ggplot(data = loc_summ, aes(x=method, y=prev)) +
  geom_boxplot() +
  ylab('Median prevalence') +
  theme_classic()

s1/s2/s3 +
  plot_annotation(tag_levels = 'a')


# species-specific prevalence
p_sum <- new_locs %>% 
  na.omit() %>% 
  group_by(method, community, species, prevalence) %>% 
  tally()

library(ggridges)
ggplot(p_sum, aes(x=prevalence, y = method)) +
  geom_density_ridges(rel_min_height=.01) +
  # xlim(0,0.6) + 
  theme_classic()

ggplot(p_sum, aes(y=prevalence, x = method, fill = method)) +
  geom_jitter(alpha = 0.2) +
  geom_violin(alpha = 0.75) +
  theme_classic() +
  coord_flip()




#### What effect does AS have on number of new locations?

new_locs %>% 
  na.omit() %>% 
  group_by(method, community) %>% 
  summarise(new_locs = length(unique(id))) %>% 
  ggplot(aes(x=method, y=new_locs)) + 
  geom_boxplot()

## that's not good - actually get more new locations using model-based approaches... Is it a problem?





#### new locations vs old locations
sp_new_locs <- new_locs %>% 
  na.omit() %>% 
  group_by(method, community, species, prevalence) %>% 
  summarise(n = sum(Observed)) %>% 
  mutate(id = paste(community, species, prevalence, sep = '_'))


sp_init_locs <- init %>% 
  na.omit() %>% 
  group_by(method, community, species, prevalence) %>% 
  summarise(n = sum(Observed)) %>% 
  mutate(id = paste(community, species, prevalence, sep = '_'))


sp_new_locs <- sp_new_locs %>% 
  rowwise() %>% 
  mutate(n_init_obs = sp_init_locs$n[match(id, sp_init_locs$id)])


ggplot(sp_new_locs, aes(x = n_init_obs, y = n)) +
  geom_point()





