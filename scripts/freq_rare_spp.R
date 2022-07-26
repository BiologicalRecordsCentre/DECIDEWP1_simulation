## working out whether rare species are more likely to be seen with 
## certain AS methods
library(tidyverse)

method = c("none", "uncertainty", "prevalence", "unc_plus_prev", "unc_plus_recs", "coverage")

AS_version = 'asv1'
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

write.csv(comm_tally, file = paste0('Outputs/', comm_version, 'communities_1km/all_observations.csv'))


obs <- read.csv('outputs/v4Community/n_observations_by_group.csv')
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




obsdf <- data.frame()
for(i in 1:50){
  print(i)
  if(file.size(fls_ls[i])==3) next() ## skip species when no observations - not ideal... 
  o<-read.csv(fls_ls[i])
  obsdf <- rbind(obsdf, o)
}

obsdf %>% 
  na.omit() %>% 
  group_by(method) %>% 
  summarise(n_spp = length(unique(species)))

obsdf %>% 
  na.omit() %>% 
  group_by(method) %>% 
  summarise(median(prevalence))

obsdf %>% 
  filter(method!='initial') %>% 
  group_by(method) %>% 
  summarise(n_seen = sum(Observed, na.rm = T))



df <- read_csv('outputs/v4Community/all_observations.csv')
head(df)

init <- df[df$method == 'initial',]

meths <- df[df$method != 'initial',]

new_locs <- meths[!meths$id %in% init$id,]

loc_summ <- new_locs %>% 
  na.omit() %>% 
  group_by(method, community, species, prevalence) %>% 
  tally() %>% 
  group_by(method, community) %>% 
  summarise(total_obs = sum(n),
            diversity = length(unique(species)))


ggplot(data = loc_summ, aes(x=method, y=total_obs)) +
  geom_boxplot()
ggplot(data = loc_summ, aes(x=method, y=diversity)) +
  geom_boxplot()




