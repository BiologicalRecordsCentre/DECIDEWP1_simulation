eval_out <- read.csv("Outputs/community_1_50_sim/community_1_50_sim_evaluation_table.csv")

#calculate differences from initial model for each method

species_list <- unique(eval_out$species)

for (i in 1:length(species_list)){
  if(!length(eval_out$mse[eval_out$species == species_list[i] & eval_out$method == "initial"]) == 0){
  eval_out$mse_diff[eval_out$species == species_list[i]] <- eval_out$mse[eval_out$species == species_list[i]]-eval_out$mse[eval_out$species == species_list[i] & eval_out$method == "initial"]
  eval_out$corr_diff[eval_out$species == species_list[i]] <- eval_out$corr[eval_out$species == species_list[i]]-eval_out$corr[eval_out$species == species_list[i] & eval_out$method == "initial"]
  eval_out$auc_diff[eval_out$species == species_list[i]] <- eval_out$auc[eval_out$species == species_list[i]]-eval_out$auc[eval_out$species == species_list[i] & eval_out$method == "initial"]
  eval_out$auc_initial[eval_out$species == species_list[i]] <- eval_out$auc[eval_out$species == species_list[i] & eval_out$method == "initial"]
  
  } else {eval_out$mse_diff[eval_out$species == species_list[i]] <- NA
  eval_out$corr_diff[eval_out$species == species_list[i]] <- NA
  eval_out$auc_diff[eval_out$species == species_list[i]] <- NA
  eval_out$auc_initial[eval_out$species == species_list[i]] <- NA}
  
}

#boxplots by method
boxplot(eval_out$mse_diff ~ eval_out$method)
boxplot(eval_out$corr_diff ~ eval_out$method)
boxplot(eval_out$auc_diff ~ eval_out$method)


#relationships with prevalence and initial model AUC
par(mar = c(5,4,4,7))
plot(eval_out$mse_diff ~ eval_out$prevalence, pch = 20, col = as.numeric(factor(eval_out$method)))
par(xpd = TRUE)
legend(x = 1.1, y = 0.1,legend = unique(eval_out$method), pch = 20, col = as.numeric(factor(eval_out$method)))
plot(eval_out$corr_diff ~ eval_out$prevalence, pch = 20, col = as.numeric(factor(eval_out$method)))
legend(x = 1.1, y = 0.1,legend = unique(eval_out$method), pch = 20, col = as.numeric(factor(eval_out$method)))
plot(eval_out$auc_diff ~ eval_out$prevalence, pch = 20, col = as.numeric(factor(eval_out$method)))
legend(x = 1.1, y = 0.1,legend = unique(eval_out$method), pch = 20, col = as.numeric(factor(eval_out$method)))

plot(eval_out$mse_diff ~ eval_out$auc_initial)
plot(eval_out$corr_diff ~ eval_out$auc_initial)
plot(eval_out$auc_diff ~ eval_out$auc_initial)

#average across species

tapply(eval_out$mse_diff, eval_out$method, function(x) mean(x, na.rm = T))
tapply(eval_out$corr_diff, eval_out$method, function(x) mean(x, na.rm = T))
tapply(eval_out$auc_diff, eval_out$method, function(x) mean(x, na.rm = T))

       