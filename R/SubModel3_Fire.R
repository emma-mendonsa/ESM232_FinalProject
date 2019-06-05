#' Add in params and other info
#'
#'

burnt_biomass_model = function(pine_biomass, fire_prob, fire_sev, 
                               forest_carbon, time = 1, carbon_coeff = 0.47) {
  
  burnt_biomass = rep(0, times=time)
  for (i in 1:time) {
    burnt_biomass = 
      ifelse(fire_prob[i] == "H" && fire_sev[i] == "Low", (sum(pine_biomass[1:3,2])*.3),
      ifelse(fire_prob[i] == "H" && fire_sev[i] == "Mod", (sum(pine_biomass[1:4,2])*.51),
      ifelse(fire_prob[i] == "H" && fire_sev[i] == "Sev", (sum(pine_biomass[1:5,2])*.8),0)))
  }
  
  post_fire_carbon = forest_carbon - (burnt_biomass * carbon_coeff)
  
  return(list(burnt_biomass, post_fire_carbon))

}



