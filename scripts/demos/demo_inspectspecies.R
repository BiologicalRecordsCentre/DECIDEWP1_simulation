### inspect individual species outputs

#dir_sp <- "" #not sure where this goes yet...I pulled outputs into local area for testing

# select community/sp of interest

spp <- "Comm9_Sp11"

# true distrbutions

community <- readRDS(paste0(dir_sp, spp, "/v3community_9_50_sim_initial.rds"))

true_dist <- community[[11]]$true_prob_occ

plot(true_dist)

# initial locations

initial_locs <- community[[11]]$observations

points(initial_locs$lat[!is.na(initial_locs$Observed)] ~ initial_locs$lon[!is.na(initial_locs$Observed)])

#initial model predictions


models <- list.files(path = paste0(dir_sp,spp,sep = "/"), pattern = paste0("(",paste( c("rf", "gam", "lr"), sep = "", collapse = "|"),")*initial.rdata"))

models_to_read <- models

#load outputs into list and extract prediction table
model_outputs <- list()
idx <- 1
for (k in models_to_read){
  try(load(paste0(paste0(dir_sp,spp,sep = "/"),k )))
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

plot(pred1)$z
