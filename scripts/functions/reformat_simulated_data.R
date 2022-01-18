

## function to reformat data into that used by cpa() function
# takes output of simulate_species() function,
# a dummy year
# and the label you want to give your species names


reformat_data <- function(sim_species_out, year=2015, species_name='Sp'){
  
  return(do.call(rbind, 
                 lapply(1:length(sim_species_out), FUN = function(i) {
                   
                   sim_species_out[[i]]$observations %>% 
                     filter(Observed==1) %>% 
                     mutate(year = year,
                            species = paste0(species_name, i))
                   
                 })))
  
}

