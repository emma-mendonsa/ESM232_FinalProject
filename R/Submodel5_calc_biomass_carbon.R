# Biomass and Carbon

calc_biomass_carbon = function(avg_stage_dbh, pop_structure, carbon_coeff = 0.47){
  
  forest_carbon_df = data.frame(avg_dbh = avg_stage_dbh,pop_structure = pop_structure)
  
  # biomass eq from Jenkins et al.
  for(i in 1:nrow(forest_carbon_df)){ # For each stage class
    
      forest_carbon_df$biomass[i] = forest_carbon_df$pop_structure[i]*exp(-2.5356 + 2.2595*log(forest_carbon_df$avg_dbh[i]))
  
    }
  
  biomass = sum(forest_carbon_df$biomass)
  carbon = biomass*carbon_coeff
  
  return(list(total_biomass=total_biomass,total_carbon=carbon))
  
}