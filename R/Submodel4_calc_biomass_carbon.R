#' Calculat the total above ground biomass and carbon in a forest
#' 
#' @param avg_stage_dbh A vector the mean or median dbh of trees in each stage
#' @param pop_structure A vector of the tree population structure
#' @param carbon_coeff The above-ground-biomass:carbon ratio 
#' @return A list of the total biomass(kg) and total carbon in a forest(kg)
#' 
#' @references Karlik and Chojnacky (2014) Biomass and Carbon data from blue oaks in a California oak savanna. Biomass and Bioenergy. 62. 228-232.

calc_biomass_carbon = function(avg_stage_dbh, pop_structure, carbon_coeff = 0.47){
  
  # Start dataframe to store results in
  forest_carbon_df = data.frame(avg_dbh = avg_stage_dbh,pop_structure = pop_structure)
  
  for(i in 1:nrow(forest_carbon_df)){ # Calculate above ground biomass in each stage class
    
      # Biomass equation published by Karlik and Chojnacky (2014)
      forest_carbon_df$biomass[i] = forest_carbon_df$pop_structure[i]*exp(-2.5356 + 2.2595*log(forest_carbon_df$avg_dbh[i]))
  
    }
  
  # Find total biomass (kg)
  biomass = sum(forest_carbon_df$biomass)
  
  # Find total carbon (kg)
  carbon = biomass*carbon_coeff
  
  return(list(total_biomass=biomass,total_carbon=carbon))
  
}