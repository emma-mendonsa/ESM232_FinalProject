#' Fire probability and severity model
#' 
#' @param initial_pop pine population structure before a fire
#' @param fire_prob probability of a fire
#' @param fire_sev severity of a fire if it occurs
#' @param forest_carbon carbon in a stand of trees/forest before a fire
#' @param low_sev Coefficient representing the effect of a low severity fire on a stand of trees/forest
#' @param mod_sev Coefficient representing the effect of a moderate severity fire on a stand of trees/forest
#' @param sev_sev Coefficient representing the effect of a high severity fire on a stand of trees/forest
#' @param time Timesteps being calculated
#' @param carbon_coeff biomass:carbon coefficient
#' 
#' @return list of: population structure after the fire, amount of biomass burnt, and carbon contained in the forest after the fire
#'

fire_model = function(initial_pop, pine_biomass, fire_prob, fire_sev, forest_carbon, low_sev = .8, mod_sev = .49, sev_sev = .2, time = 1, carbon_coeff = 0.47) {
  
  postFire_pop = ifelse(fire_prob =="L", list(c(initial_pop[1], # If the probability of a fire is low, we assume there is no fire. 
                                 initial_pop[2],
                                 initial_pop[3],
                                 initial_pop[4],
                                 initial_pop[5])),
             
               ifelse(fire_prob == "H" && fire_sev == "Low", # High probability of fire, and the effect of low fire severity
                     list(c(initial_pop[1]*low_sev,
                       initial_pop[2]*low_sev,
                       initial_pop[3]*low_sev,
                       initial_pop[4],
                       initial_pop[5])),
                     
                     ifelse(fire_prob == "H" && fire_sev == "Mod", # High probability of fire, and the effect of moderate fire severity
                            list(c(initial_pop[1]*mod_sev,
                              initial_pop[2]*mod_sev,
                              initial_pop[3]*mod_sev,
                              initial_pop[4]*mod_sev,
                              initial_pop[5])),
                            
                            ifelse(fire_prob == "H" && fire_sev == "Sev" && pine_biomass > 10000, # High probability of fire, and the effect of severe fire
                                   list(c(initial_pop[1]*sev_sev,
                                     initial_pop[2]*sev_sev,
                                     initial_pop[3]*sev_sev,
                                     initial_pop[4]*sev_sev,
                                     initial_pop[5]*sev_sev)),
                                   
                                   return("Unknown fire model input"))))) # If some other combination return an error
   
  # Create empty data.frame for change in biomass output                                                                
  burnt_biomass = rep(0, times=time) 
  
  for (i in 1:time) {
    burnt_biomass = 
      ifelse(fire_prob[i] == "H" && fire_sev[i] == "Low", (pine_biomass*.2), # If there is a low severity firt, 20% of biomass will burn
      ifelse(fire_prob[i] == "H" && fire_sev[i] == "Mod", (pine_biomass*.51), # If there is a moderate severity firt, 51% of biomass will burn
      ifelse(fire_prob[i] == "H" && fire_sev[i] == "Sev", (pine_biomass*.8),0))) # If there is a high severity firt, 80% of biomass will burn
  }
  
  post_fire_carbon = forest_carbon - (burnt_biomass * carbon_coeff) # Calculate the carbon in the forest after the fire
  
  return(list(postFire_pop = unlist(postFire_pop), # Population structure after the fire
              burnt_biomass = burnt_biomass, # Burnt biomass (kg)
              post_fire_carbon = post_fire_carbon)) # Carbon remaining in the forest after the fire (kg)

}



