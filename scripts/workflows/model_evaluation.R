eval_out <- read.csv("Outputs/community_4_50_sim/community_4_50_sim_evaluation_table.csv")

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


png("Community 4 results by method.png", height = 800, width = 400)
par(mfrow=c(3,1))
#boxplots by method
boxplot(eval_out$mse_diff ~ eval_out$method, xlab = "Method", ylab = "MSE difference from initial model")
boxplot(eval_out$corr_diff ~ eval_out$method, xlab = "Method", ylab = "Correlation difference from initial model")
boxplot(eval_out$auc_diff ~ eval_out$method, xlab = "Method", ylab = "AUC difference from initial model")
dev.off()

#relationships with prevalence and initial model AUC

png("Community 4 results against prevalence.png", height = 1600, width = 1200, res = 200, pointsize = 16)
par(mfrow=c(3,1))
par(mar = c(4,4,1,7))
plot(eval_out$mse_diff ~ eval_out$prevalence, pch = 20, col = as.numeric(factor(eval_out$method)), xlab = "Prevalence", ylab = "MSE difference from initial model")
par(xpd = TRUE)
legend(x = 1.1, y = 0.1,legend = unique(eval_out$method), pch = 20, col = as.numeric(factor(eval_out$method)))
plot(eval_out$corr_diff ~ eval_out$prevalence, pch = 20, col = as.numeric(factor(eval_out$method)), xlab = "Prevalence", ylab = "Correlation difference from initial model")
legend(x = 1.1, y = 0.1,legend = unique(eval_out$method), pch = 20, col = as.numeric(factor(eval_out$method)))
plot(eval_out$auc_diff ~ eval_out$prevalence, pch = 20, col = as.numeric(factor(eval_out$method)), xlab = "Prevalence", ylab = "AUC difference from initial model")
legend(x = 1.1, y = 0.1,legend = unique(eval_out$method), pch = 20, col = as.numeric(factor(eval_out$method)))
dev.off()


png("Community 4 results against AUC.png", height = 1600, width = 1200, res = 200, pointsize = 16)
par(mfrow=c(3,1))
par(mar = c(4,4,1,7))
plot(eval_out$mse_diff ~ eval_out$auc_initial, pch = 20, col = as.numeric(factor(eval_out$method)), xlab = "AUC of initial model", ylab = "MSE difference from initial model")
par(xpd = TRUE)
legend(x = 1.1, y = 0.1,legend = unique(eval_out$method), pch = 20, col = as.numeric(factor(eval_out$method)))
plot(eval_out$corr_diff ~ eval_out$auc_initial, pch = 20, col = as.numeric(factor(eval_out$method)), xlab = "AUC of initial model", ylab = "Correlation difference from initial model")
legend(x = 1.1, y = 0.1,legend = unique(eval_out$method), pch = 20, col = as.numeric(factor(eval_out$method)))
plot(eval_out$auc_diff ~ eval_out$auc_initial, pch = 20, col = as.numeric(factor(eval_out$method)), xlab = "AUC of initial model", ylab = "AUC difference from initial model")
legend(x = 1.1, y = 0.1,legend = unique(eval_out$method), pch = 20, col = as.numeric(factor(eval_out$method)))
dev.off()




#average across species

tapply(eval_out$mse_diff, eval_out$method, function(x) mean(x, na.rm = T))
tapply(eval_out$corr_diff, eval_out$method, function(x) mean(x, na.rm = T))
tapply(eval_out$auc_diff, eval_out$method, function(x) mean(x, na.rm = T))


library(tidyverse)
library(ggridges)

# outputs from susan’s first run.
cdf <- eval_out
head(cdf)


ggplot(cdf, aes(x = prevalence, y = corr, colour = method)) +
  geom_point() +
  geom_smooth()


cdf2 <- cdf %>% 
  pivot_longer(cols = 3:5)


ggplot(cdf2, aes(x = prevalence, y = value, colour = method)) +
  geom_point() +
  geom_smooth(se = F) +
  facet_wrap(~name) +
  labs(y='Evaluation metric') +
  theme_bw()

# ggplot(cdf2, aes(x = value, y = method, fill = method)) +
#   geom_density_ridges()+
#   facet_wrap(~name) +
#   labs(y='Evaluation metric') +
#   theme_bw()

ggplot(cdf2, aes(x = method, y = value, colour = name)) +
  geom_boxplot() +
  facet_wrap(~name, scales = 'free') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

       