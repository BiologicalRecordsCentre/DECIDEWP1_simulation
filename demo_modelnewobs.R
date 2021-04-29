#' # Run model with new adaptively sampled observations and compare to previous model
#' 
#' 
library(raster)
library(virtualspecies)
library(dismo)
library(rgdal)
library(knitr)

set.seed(1000)

#load output of demo_simulatenewobs.R 
load("virt_comm_10spp_AS_newdata.Rdata")

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

hbv_y <- raster::stack("hbv_y.grd") 

#for running the models we subset the environmental data to the variables used in the species generation
env_data <- subset(hbv_y, subset = virt_comm1[[1]]$details$variables)

#source code from Thomas' workflow
source("getpredictions.R")

#' Example for Sp1 (not on LOTUS, having issues today)

species = "Sp1"
model <- "lr"
k = 4

#run model for first species
sdm <- fsdm(species = species, model = "lr",
            climDat = env_data, spData = pa_sets, knots_gam = -1,
            k = 4, 
            write =  TRUE, outPath = "lr_outs/AS_newdata/")

#predictions

preds1 <- get_predictions(sdm, "lr", env_data)


## save files ##
species_name <- gsub(pattern = ' ', replacement = '_', species) # get species name without space

outPath <- paste0(getwd(), "/Outputs/AS_newdata/")

# save prediction raster
writeRaster(x = preds1$mean_predictions, 
            filename = paste0(outPath, model, "_SDMs_", species_name, "_meanpred.grd"),
            format = 'raster', overwrite = T)

# save quantile max min
writeRaster(x = preds1$quant_minmax, 
            filename = paste0(outPath, model, "_SDMs_", species_name, "_quantilemaxmin.grd"),
            format = 'raster', overwrite = T)

# save quantile range raster
writeRaster(x = preds1$quant_range, 
            filename = paste0(outPath, model, "_SDMs_", species_name, "_quantilerange.grd"),
            format = 'raster', overwrite = T)

# write AUC to file for easy-access
write.csv(x = data.frame(raw_AUC = sdm$AUC,
                         meanAUC = sdm$meanAUC),
          file = paste0(outPath, model, "_SDMs_", species_name, "_AUC_values.csv"))

# write data to file too
write.csv(x = sdm$Data,
          file = paste0(outPath, model, "_SDMs_", species_name, "_Data.csv"))

# remove data from model output
sdm$Data <- NULL

# output of model to store
model_output <- list(species = species_name,
                     model = model,
                     sdm_output = sdm,
                     number_validations = k)

save(model_output, file = paste0(outPath, model, "_SDMs_", species_name, 
                                 ".rdata"))



#' Plot new model outputs alongside original outputs
#' 

#load original model predictions and quantile range
orig_preds <- raster("Outputs/lr_SDMs_Sp1_meanpred.grd")
orig_rnge <- raster("Outputs/lr_SDMs_Sp1_quantilerange.grd")
orig_score <- raster("Outputs/lr_SDMs_Sp1_DECIDEscore.grd")

png("Plots/Sp1_AS_example.png", height = 200, width = 300, res = 300, units = "mm", pointsize = 12)

par(mfrow=c(3,3))
par(mar = c(2,2,2,2))
plot(virt_comm1[[1]], main = "Environmental suitability")
plot(pa[[1]]$probability.of.occurrence, main = "Probability of occurrence")
plot(pa[[1]]$pa.raster, main = "Presence absence")
points(sp.obs[[1]]$sample.points[is.na(sp.obs[[1]]$sample.points$Observed),1:2], pch = 20)
plot(orig_score, main = "DECIDE score - focal species")
points(new.obs[[1]]$sample.points[is.na(new.obs[[1]]$sample.points$Observed),1:2], pch = 20)
plot(DECIDE_allspp, main = "DECIDE score - all species")
points(new.obs[[1]]$sample.points[is.na(new.obs[[1]]$sample.points$Observed),1:2], pch = 20)
plot(orig_preds, main = "Original predicted prob. occ")
plot(orig_rnge, main = "Original quantile range of predictions")
plot(preds1$mean_predictions, main = "Predicted prob. occ after adaptive sampling")
plot(preds1$quant_range, main = "Quantile range of predictions after adaptive sampling")


dev.off()

#' ## Evaluate models
#' 
#' We can calculate the RMSE and correlation between the true probability of occurrence and the estimated probabilities for the original model (before adaptive sampling) and the revised model (after adaptive sampling)
#' 

rmse_orig <- sqrt(mean(getValues((pa[[1]]$probability.of.occurrence-orig_preds)^2), na.rm=T))

rmse_AS <- sqrt(mean(getValues((pa[[1]]$probability.of.occurrence-preds1$mean_predictions)^2), na.rm=T))

rmse_orig
rmse_AS
#new data has slightly reduced RMSE compared to original model

corr_orig <- cor(getValues(pa[[1]]$probability.of.occurrence), getValues(orig_preds),use = "complete.obs")#may or may not be calculating this correctly

corr_AS <- cor(getValues(pa[[1]]$probability.of.occurrence), getValues(preds1$mean_predictions),use = "complete.obs")#may or may not be calculating this correctly

corr_orig
corr_AS
#new data slightly increases correlation between true and estimated probabilities