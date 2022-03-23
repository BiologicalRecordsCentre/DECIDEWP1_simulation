### inspect individual species outputs

#dir_sp <- "" #not sure where this goes yet...I pulled outputs into local area for testing

# select community/sp of interest

comm <- "11"

spp <- "24"

#set up to compare two different AS methods side by side

method1 <- "none"

method2 <- "unc_plus_recs"

# true distrbutions

community <- readRDS(paste0(dir_sp, "Comm",comm, "_Sp",spp, "/v3community_",comm,"_50_sim_initial.rds"))

true_dist <- community[[as.numeric(spp)]]$true_prob_occ

plot(true_dist)

# initial locations

initial_locs <- community[[as.numeric(spp)]]$observations

points(initial_locs$lat[!is.na(initial_locs$Observed)] ~ initial_locs$lon[!is.na(initial_locs$Observed)])

#initial model predictions


models <- list.files(path = paste0(dir_sp,"Comm",comm, "_Sp",spp,"/"), pattern = paste0("(",paste( c("rf", "gam", "lr"), sep = "", collapse = "|"),")*initial.rdata"))

pred_func <- function(models){
  
models_to_read <- models

#load outputs into list and extract prediction table
model_outputs <- list()
idx <- 1
for (k in models_to_read){
  try(load(paste0(paste0(dir_sp,"Comm",comm, "_Sp",spp,"/"),k )))
  model_type <- model_output$model
  if (model_type != "rf"){
    model_preds <- model_output$predictions}
  if (model_type == "rf"){
    model_preds <- model_output$predictions[,names(model_output$predictions) %in% c("x", "y", "mean.1", "sd", "DECIDE_score.1")]
    names(model_preds) <- c("x", "y", "mean", "sd", "DECIDE_score")
  }
  model_preds$mean <- model_preds$mean#*(1-prevalence_vec[j]) #weight prevalence by rarity
  model_preds$DECIDE_score <- model_preds$mean*model_preds$sd #recalculate with prevalence weighted by rarity
  model_outputs[[idx]] <- model_preds
  names(model_outputs)[idx] <- model_type
  idx <- idx + 1
}

#average model outputs (note - not weighted by AUC currently)
mod_average <- Reduce(`+`, model_outputs) / length(model_outputs)

names(mod_average)[3] <- "z"

pred1 <- rasterFromXYZ(mod_average)

return(pred1)

}

pred1 <- pred_func(models)

par(mfrow=c(1,2))
plot(pred1$z)#initial predictions
plot(pred1$sd)#initial uncertainty

#method 1

AS1 <- readRDS(paste0(dir_sp, "Comm",comm, "_Sp",spp, "/asv8_v3community_",comm,"_50_sim_initial_AS_",method1,".rds"))

AS_locs1 <- AS1[[as.numeric(spp)]]$observations

plot(true_dist, main = method1)
points(initial_locs$lat[!is.na(initial_locs$Observed)] ~ initial_locs$lon[!is.na(initial_locs$Observed)], pch = 20)
points(AS_locs1$lat[!is.na(AS_locs1$Observed)] ~ AS_locs1$lon[!is.na(AS_locs1$Observed)], col = "red")

#method 2

AS2 <- readRDS(paste0(dir_sp, "Comm",comm, "_Sp",spp, "/asv8_v3community_",comm,"_50_sim_initial_AS_",method2,".rds"))

AS_locs2 <- AS2[[as.numeric(spp)]]$observations

plot(true_dist, main = method2)
points(initial_locs$lat[!is.na(initial_locs$Observed)] ~ initial_locs$lon[!is.na(initial_locs$Observed)], pch = 20)
points(AS_locs2$lat[!is.na(AS_locs2$Observed)] ~ AS_locs2$lon[!is.na(AS_locs2$Observed)], col = "red")

#coverage has no new locations sampled...data going into models is the same, eval differences just a result of stochasiticity in cross validation

#predictions

#method 1
modelsAS1 <- list.files(path = paste0(dir_sp,"Comm",comm, "_Sp",spp,"/"), pattern = paste0("(",paste( c("rf", "gam", "lr"), sep = "", collapse = "|"),")*initial_AS_",method1,".rdata"))

predAS1 <- pred_func(modelsAS1)

par(mfrow=c(1,2))
plot(predAS1$z)#initial predictions
plot(predAS1$sd)#initial uncertainty

#method 2
modelsAS2 <- list.files(path = paste0(dir_sp,"Comm",comm, "_Sp",spp,"/"), pattern = paste0("(",paste( c("rf", "gam", "lr"), sep = "", collapse = "|"),")*initial_AS_",method2,".rdata"))

predAS2 <- pred_func(modelsAS2)

par(mfrow=c(1,2))
plot(predAS2$z)#initial predictions
plot(predAS2$sd)#initial uncertainty

#compare predicted maps side by side

par(mfrow=c(2,2))
plot(true_dist, main = "true distribution")
plot(pred1$z, main = "initial")
plot(predAS1$z, main = method1)
plot(predAS2$z, main = method2)
