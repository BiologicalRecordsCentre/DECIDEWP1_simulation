err.files <- list.files("_rslurm_sdm_simulated_species/", pattern = ".err")

pars <- read.csv("_rslurm_sdm_simulated_species/pars.csv")

err_out <- list()

for (i in 1:length(err.files)){
  error_log <- readLines(paste0("_rslurm_sdm_simulated_species/",err.files[i]), skipNul = T)
  err_out[[i]] <- error_log[grep(pattern= "error|Error", error_log)]

}

pars$error_msg <- unlist(lapply(err_out,function(x) if(identical(x,character(0))) ' ' else x))

batchid <- unique(pars$BatchID)

write.csv(pars, paste0(getwd(), "/Outputs/logs/log_", batchid,".csv"))