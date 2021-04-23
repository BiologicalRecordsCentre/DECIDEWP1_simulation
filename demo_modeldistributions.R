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
library(dismo)

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
pa_sets[[i]] <- list(pa_data = pa_data, Presence = pa_data[pa_data$Observed == 1,1:2], pseudoAbsence = pa_data[pa_data$Observed == 0,1:2])
names(pa_sets)[i] <- paste0("Sp",i)
}

#' Pseudoabsence weighting dependent on model, Thomas' code has this accounted for


#' 3. Run DECIDE models 
#' 
#' Initially we can run just the logistic regression models as a test

#source Edited Rob Functions
devtools::source_url("https://raw.githubusercontent.com/TMondain/DECIDE_WP1/main/scripts/functions/Edited_Rob_Functions.R")

#run model for first species
lr_mod <- fsdm(species = "Sp1", model = "lr",
     climDat = subset(hbv_y, subset = virt_comm1[[1]]$details$variables), spData = pa_sets, knots_gam = -1,
     k = 4, 
     write =  TRUE, outPath = "lr_outs/")

lr_mod




