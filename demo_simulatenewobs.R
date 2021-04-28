#' # Combine DECIDE scores across species and derive "nudges"
#' 
#' We can use methods derived for the real DECIDE sampling to do this
#' 
#' Firstly, read in the species model outputs - we need the DECIDE score for each species (skipping the decision step on this for now)
#' 
#' 

library(raster)
library(virtualspecies)

set.seed(1000)

filenames <- list.files("Outputs", pattern = "*DECIDEscore.grd", full.names = TRUE)

scores <- raster::stack(filenames)

#' Read in functions from DECIDE_adaptive sampling repo
#' 
devtools::source_url("https://raw.githubusercontent.com/BiologicalRecordsCentre/DECIDE_adaptivesampling/main/Scripts/modules/recommend_rank.R")

# doesn't work 
#r1 <- recommend_rank(scores, method = "additive")

#' Skip ranking step for now, just get the average DECIDE score across species as placeholder for aggegrating across species
#' 
DECIDE_allspp <- calc(scores, fun = mean, na.rm=T)
par(mfrow=c(1,1))
plot(DECIDE_allspp)

#' Use this layer to determine new sampling bias
#' 
# load simulated species distributions

load("virt_comm_10spp.Rdata")

#sample new obs (no longer weighted by prevalence? want to weight by prevalence? not sure) weighted by DECIDE score

DECIDE_weights <- DECIDE_allspp*100 #rescale (needs refining)

par(mfrow=c(5,2))
new.obs <- list()
for (i in 1:10){
  max_obs <- round(prev_vec[i]*100)
  new.obs[[i]] <- sampleOccurrences(pa[[i]], n = max_obs, type = "presence only", detection.probability = 0.5, bias = "manual", weights = DECIDE_weights)
  names(new.obs[[i]]$sample.points) <- c("lon", "lat", "Real", "Observed")
}



