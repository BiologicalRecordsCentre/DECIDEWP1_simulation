### inspect individual species outputs

#dir_sp <- "" #not sure where this goes yet...I pulled outputs into local area for testing

# select community/sp of interest

comm <- "9"

spp <- "11"

method1 <- "none"

method2 <- "coverage"

# true distrbutions

community <- readRDS(paste0(dir_sp, "Comm",comm, "_Sp",spp, "/v3community_",comm,"_50_sim_initial.rds"))

true_dist <- community[[11]]$true_prob_occ

plot(true_dist)

# initial locations

initial_locs <- community[[11]]$observations

points(initial_locs$lat[!is.na(initial_locs$Observed)] ~ initial_locs$lon[!is.na(initial_locs$Observed)])

#initial model predictions


models <- list.files(path = paste0(dir_sp,"Comm",comm, "_Sp",spp,"/"), pattern = paste0("(",paste( c("rf", "gam", "lr"), sep = "", collapse = "|"),")*initial.rdata"))

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


