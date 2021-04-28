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

DECIDE_weights <- DECIDE_allspp

#scale to probability raster - sum of all cells = 1 so probs are small...could cause computational issues?
DECIDE_weights <- DECIDE_weights/sum(getValues(DECIDE_weights), na.rm=T)

#sample new observations with the same probability for all species

#function from https://scrogster.wordpress.com/2012/04/22/using-r-for-spatial-sampling-with-selection-probabilities-defined-in-a-raster/
probsel<-function(probrast, N){
  x<-getValues(probrast)
  #set NA cells in raster to zero
  x[is.na(x)]<-0
  samp<-sample(nrow(probrast)*ncol(probrast), size=N, prob=x)
  samprast<-raster(probrast)
  samprast[samp]<-1 #set value of sampled squares to 1
  #convert to SpatialPoints
  points<-rasterToPoints(samprast, fun=function(x){x>0})
  points<-SpatialPoints(points)
  return(points)
}

#new sampling locations
new_locs <- probsel(DECIDE_weights, 100)

#now sample from species rasters with given locations
new.obs <- list()
for (i in 1:length(pa)){
  sample.points <- data.frame(coordinates(new_locs)[,1:2])
  sample.points$Real <- extract(pa[[i]]$pa.raster, sample.points)
  sample.points <- sample.points[sample.points$Real == 1,]#remove absences to create presence-only data?
  sample.points$Observed <- sample.points$Real * (rbinom(nrow(sample.points),1,0.5))
  sample.points[sample.points == 0] <- NA
  new.obs[[i]] <- list()
  new.obs[[i]]$sample.points <- sample.points
}

#should mimic structure of sampleOccurrences output to help with modelling although some weird NA locations need checking, also may be simpler way to do this

