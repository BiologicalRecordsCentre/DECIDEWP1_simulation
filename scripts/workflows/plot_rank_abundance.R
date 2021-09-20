community1 <- readRDS("results_0.RDS")[[1]]
community2 <- readRDS("results_1.RDS")[[1]]

rank_abundance <- function(community, add = FALSE){

prev_vec <- vector()
for (i in 1:length(community)){
  prev_vec[i] <- cellStats(community[[i]]$pres_abs, sum)/ncell(community[[i]]$pres_abs)
}
prev_vec

rank_prev <- sort(prev_vec, decreasing = TRUE)

if(add == FALSE){
plot(rank_prev, type= "l", xlab = "Abundance rank", ylab = "Prevalence")}

if(add == TRUE){
  lines(rank_prev, col = sample(1:10,1))
}
}

rank_abundance(community1)
rank_abundance(community2, add = TRUE)