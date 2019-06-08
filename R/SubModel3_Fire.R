#' Add in params and other info
#'
#'

burnt_biomass_model = function(initial_pop, pine_biomass, fire_prob, fire_sev, forest_carbon, low_sev = .7, mod_sev = .49, sev_sev = .2, time = 1, carbon_coeff = 0.47) {
  
  
  postFire_pop = ifelse(fire_prob =="L", list(c(initial_pop[1],
                                 initial_pop[2],
                                 initial_pop[3],
                                 initial_pop[4],
                                 initial_pop[5])),
             
               ifelse(fire_prob == "H" && fire_sev == "Low", 
                     list(c(initial_pop[1]*low_sev,
                       initial_pop[2]*low_sev,
                       initial_pop[3]*low_sev,
                       initial_pop[4],
                       initial_pop[5])),
                     
                     ifelse(fire_prob == "H" && fire_sev == "Mod", 
                            list(c(initial_pop[1]*mod_sev,
                              initial_pop[2]*mod_sev,
                              initial_pop[3]*mod_sev,
                              initial_pop[4]*mod_sev,
                              initial_pop[5])),
                            
                            ifelse(fire_prob == "H" && fire_sev == "Sev", 
                                   list(c(initial_pop[1]*sev_sev,
                                     initial_pop[2]*sev_sev,
                                     initial_pop[3]*sev_sev,
                                     initial_pop[4]*sev_sev,
                                     initial_pop[5]*sev_sev)),NA))))
  
  
  burnt_biomass = rep(0, times=time)
  
  for (i in 1:time) {
    burnt_biomass = 
      ifelse(fire_prob[i] == "H" && fire_sev[i] == "Low", (pine_biomass*.3),
      ifelse(fire_prob[i] == "H" && fire_sev[i] == "Mod", (pine_biomass*.51),
      ifelse(fire_prob[i] == "H" && fire_sev[i] == "Sev", (pine_biomass*.8),0)))
  }
  
  post_fire_carbon = forest_carbon - (burnt_biomass * carbon_coeff)
  
  return(list(postFire_pop = unlist(postFire_pop), burnt_biomass = burnt_biomass, post_fire_carbon = post_fire_carbon))

}



