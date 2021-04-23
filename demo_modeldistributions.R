#' # Example of basic modelling of simulated species using northern England region
#' 
#' The steps shown in this document are as follows:
#' 
#' 1. Simulate species distributions (following steps outlined in demo_simulatebaseline)
#' 
#' 2. Create input data for models (create pseudo-absences)
#' 
#' 3. Fit DECIDE SDMs to the distributions
#' 
#' 4. Evaluate DECIDE score across region
#' 
#' 
#' 
#' ## 1. Simulate distributions (or read in simulated spp)
library(raster)
library(virtualspecies)

set.seed(1000)

#load output of demo_simulatebaseline.R (could be reduced in size to speed this up)
load("virt_comm_10spp.Rdata")

#' ## 2. Create input data for models
#' 
#' We need to extract the virtual species data we simulated and combine into a community dataset. For each species we can create a new dataset with pseudoabsences generated from the other species in the community
#' 
#create a pseudo-absence dataset
pa_sets <- list()

for(i in 1:length(sp.obs)){
  #extract occurrence records
occs <- sp.obs[[i]]$sample.points[!is.na(sp.obs[[i]]$sample.points$Observed),]
abs <- data.frame()
  for (j in 1:length(sp.obs)){
    if(j != i) {
      #observations of other species treated as absences for species of interest
      abs <- rbind(abs, sp.obs[[j]]$sample.points[!is.na(sp.obs[[j]]$sample.points$Observed),])
    }
  }
abs$Observed <- 0
#check no duplicates
pa_data <- rbind(occs, abs)
dups <- pa_data[duplicated(pa_data[1:2]),]
#always remove duplicates - because occurrences bound first then any absences in same location as occurrences will be removed, also removed any duplicates in pseudo absence data
if(nrow(dups) == 0) {next} else {
  pa_data <- pa_data[!duplicated(pa_data[1:2]),]
}
pa_sets[[i]] <- pa_data
}

#' Maybe need to consider pseudo-absence weighting...to come back to 
#' 
#' 



