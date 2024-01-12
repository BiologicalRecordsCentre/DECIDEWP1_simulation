
library(rslurm)

dirs <- config::get("LOTUSpaths")

source('../scripts/slurm_run_sim_sdm_function.R')

# name of the versions we are running - so we're not overwriting things
# one for community-level which includes the community folders and species models folders
community_version = 'v4'


## run multiple adaptive sampling versions at once
asv_version = data.frame(as_version =  c("asv1", "asv2", "asv3", "asv4"),
                         uptake_value = c(0.1, 0.01, 0, 0.5))

for(asv in 1:nrow(asv_version)){
  
  # One name for the adaptive sampling round to allow us to create different sampling methods of the same initial community
  # This doesn't get used if running only the initial model, but does get used when running the adaptive sampling methods.
  AS_version = asv_version$as_version[asv]
  print(AS_version)
  
  # the name of the simulation run - same as slurm_simulate species
  simulation_run_name = 'communities_1km'
  
  ## create for loop to submit all scripts at once
  job_seqs <- list(1:11, 12:22, 23:33, 34:44, 45:50) 
  
  for(s in job_seqs) {
    
    ## new parameters code to try and automate the parameter generation file a little more
    n_species = 1:50 # vector of number of species in each community
    n_communities = s # number of communities to go through # can submit 11 at a time ## done 1-11, running 12-22, done 23-33, done 34-44, done 45-50
    models = c('lr', 'gam', 'rf') # models to run
    data_type = "initial_AS_coverage" #c("initial_AS_none", "initial_AS_uncertainty", "initial_AS_prevalence", "initial_AS_unc_plus_prev", "initial_AS_unc_plus_recs", "initial_AS_coverage") # 'initial'
    
    # name of the versions we are running - so we're not overwriting things
    # one for community-level which includes the community folders and species models folders
    community_version = community_version
    
    # One name for the adaptive sampling round to allow us to create different sampling methods of the same initial community
    # This doesn't get used if running only the initial model, but does get used when running the adaptive sampling methods.
    AS_version = AS_version
    
    # the name of the simulation run - same as slurm_simulate species
    simulation_run_name = simulation_run_name
    
    # # set commpath to something for testing (then delete 'dirs$')
    # dirs <- data.frame(commpath = 'blob')
    
    # pars data frame - written to be able to automatically generate the correct pars file depending on inputs above.
    pars <- data.frame(index = rep(n_species, length(n_communities)*length(models)*length(data_type)),
                       spdata = rep(sprintf(
                         paste0(dirs$commpath, community_version, simulation_run_name,"/", community_version,
                                "community_%i_%i_sim/", ifelse(data_type!='initial', paste0(AS_version, '_'), ''), community_version, "community_%i_%i_sim_%s.rds"), 
                         rep(rep(n_communities, each = length(models)), each = length(data_type)), max(n_species), 
                         rep(rep(n_communities, each = length(models)), each = length(data_type)), max(n_species), data_type
                       ), each = length(n_species)), # location of all the community data
                       model = rep(rep(rep(models, length(n_communities)), each = length(data_type)), each = length(n_species)), # which models to use
                       data_type = rep(rep(data_type, length(n_communities)*length(models)), each = length(n_species)), 
                       writeRas = FALSE,
                       GB = TRUE,
                       community_version = community_version,
                       simulation_run_name = simulation_run_name,
                       AS_version = AS_version,
                       n_communities = rep(rep(rep(n_communities, each = length(models)), each = length(data_type)), each = length(n_species)),
                       n_species = max(n_species)
    )
    
    print(dim(pars))
    
    #test with subset of runs
    #pars <- pars[-c(1,78,103,166,208,294,315,352,422, 484,517, 569, 646,675, 735),]
    
    # # resubmit long runs - job numbers from lotus
    # resub_rows <- c(1486, 1485, 1455, 1456, 1388, 1389, 1395, 1396, 1365, 1366, 1335, 1336, 1305, 1306, 1245, 1246, 1238, 1239)+1
    # pars <- pars[resub_rows,]
    
    #### slurm apply call
    sdm_slurm <- slurm_apply(slurm_run_sim_sdm,
                             params = pars,
                             jobname = paste0(community_version, '_', AS_version, '_sdm_simulated_species_', max(n_communities)),
                             nodes = length(pars$index),
                             cpus_per_node = 1,
                             slurm_options = list(partition = 'short-serial',#-4hr',
                                                  time = '23:59:59',
                                                  mem = 10000,
                                                  output = "sim_sdm_%a.out",
                                                  error = "sim_sdm_%a.err"),#,
                             # account = "short4hr"),
                             sh_template = "jasmin_submit_sh.txt",
                             submit = T)
    pars$BatchID <- sdm_slurm$jobid
    pars$JobID <- 0:(nrow(pars)-1)#slurm job ID
    write.csv(pars, paste0('_rslurm_', community_version, '_', AS_version, '_sdm_simulated_species_', max(n_communities),'/pars.csv'))#to match to error files
    
  }
  
}
