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

sdm <- lr_mod
#predictions

env_data <- subset(hbv_y, subset = virt_comm1[[1]]$details$variables)

boots_out <- raster::stack(lapply(lr_mod$Bootstrapped_models, FUN = function(x) predict(env_data, x, type="response", index=NULL)))

## quantiles
print(paste0('#####   getting quantiles   #####'))
mean_preds <- calc(boots_out, fun = mean, na.rm = T)
quant_preds <- calc(boots_out, fun = function(x) {quantile(x, probs = c(0.05, 0.95), na.rm = TRUE)})
rnge <- quant_preds[[2]]-quant_preds[[1]]

## save files ##
print("#####     Saving files     #####")
species_name <- gsub(pattern = ' ', replacement = '_', species) # get species name without space

outPath <- paste0(getwd(), "/Outputs/")

# save prediction raster
print("#####     Saving prediction raster     #####")
writeRaster(x = mean_preds, 
            filename = paste0(outPath, model, "_SDMs_", species_name, "_meanpred.grd"),
            format = 'raster', overwrite = T)

# save quantile max min
print("#####     Saving quantile max min raster     #####")
writeRaster(x = quant_preds, 
            filename = paste0(outPath, model, "_SDMs_", species_name, "_quantilemaxmin.grd"),
            format = 'raster', overwrite = T)

# save quantile range raster
print("#####     Saving quantile range raster     #####")
writeRaster(x = rnge, 
            filename = paste0(outPath, model, "_SDMs_", species_name, "_quantilerange.grd"),
            format = 'raster', overwrite = T)

# write AUC to file for easy-access
print("#####     Writing AUC to file     #####")
write.csv(x = data.frame(raw_AUC = sdm$AUC,
                         meanAUC = sdm$meanAUC),
          file = paste0(outPath, model, "_SDMs_", species_name, "_AUC_values.csv"))

# write data to file too
print("#####     Writing data to file     #####")
write.csv(x = sdm$Data,
          file = paste0(outPath, model, "_SDMs_", species_name, "_Data.csv"))

# save subset model output
print("#####     Saving model output     #####")

# remove data from model output
sdm$Data <- NULL

# outout of model to store
model_output <- list(species = species_name,
                     model = model,
                     sdm_output = sdm,
                     number_validations = k)

save(model_output, file = paste0(outPath, model, "_SDMs_", species_name, 
                                 ".rdata"))

#' Calculate very simple DECIDE score - prediction * quantile range

DECIDE_score <- mean_preds*rnge

#' Plot maps
#' 
par(mfrow=c(3,2))
par(mar = c(2,2,2,2))
plot(virt_comm1[[1]], main = "Environmental suitability")
plot(pa[[1]]$probability.of.occurrence, main = "Probability of occurrence")
plot(pa[[1]]$pa.raster, main = "Presence absence")
points(sp.obs[[1]]$sample.points[is.na(sp.obs[[1]]$sample.points$Observed),1:2], pch = 20)
plot(mean_preds, main = "Predicted prob. occ")
plot(rnge, main = "Quantile range of predictions")
plot(DECIDE_score, main = "DECIDE score")


