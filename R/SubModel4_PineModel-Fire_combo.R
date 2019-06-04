#' Add in params and other info
#'
#'
#'
#'
#'
#'
#'
#'
#'

whitePine_matrix_FIREmodel = function(fertility, survival, growth, 
                                  initial_pop, fire_prob, fire_sev, time, carbon_coeff = 0.47) {
  
  nstages = length(fertility)
  
  # Initialize Leslie Matrix
  whitePine_matrix = matrix(nrow=nstages, ncol=nstages)
  colnames(whitePine_matrix) = stages
  rownames(whitePine_matrix) = stages
  
  # Fill in with zero values
  whitePine_matrix[,] = 0.0
  
  # Put fertility vector in the first row
  whitePine_matrix[1,] = fertility
  
  # Add in survival values. This is particularly for survival & growth -- the probability that an individual survivies and grows into the next stage
  for (i in 1:(nstages)) {
    whitePine_matrix[i,i] = survival[i]
  }
 
  for (i in 1:(nstages-1)) {
    whitePine_matrix[i+1,i] = growth[i]
  } 

#Fire severity coefficients:
  low_sev = .7
  mod_sev = .49
  sev_sev = .2
  
  stable_stage = as.numeric(stable.stage(whitePine_matrix))
  P0 = ifelse(fire_prob =="L", c(initial_pop*stable_stage[1],initial_pop*stable_stage[2],
                                 initial_pop*stable_stage[3],initial_pop*stable_stage[4],
                                 initial_pop*stable_stage[5]),
              ifelse(fire_prob == "H" && fire_sev == "Low", 
                     c(initial_pop*low_sev*stable_stage[1],
                       initial_pop*low_sev*stable_stage[2],
                       initial_pop*low_sev*stable_stage[3],
                       initial_pop*stable_stage[4],
                       initial_pop*stable_stage[5]),
              ifelse(fire_prob == "H" && fire_sev == "Mod", 
                     c(initial_pop*mod_sev*stable_stage[1],
                       initial_pop*mod_sev*stable_stage[2],
                       initial_pop*mod_sev*stable_stage[3],
                       initial_pop*mod_sev*stable_stage[4],
                       initial_pop*stable_stage[5]),
              ifelse(fire_prob == "H" && fire_sev == "Sev", 
                     c(initial_pop*sev_sev*stable_stage[1],
                       initial_pop*sev_sev*stable_stage[2],
                       initial_pop*sev_sev*stable_stage[3],
                       initial_pop*sev_sev*stable_stage[4],
                       initial_pop*sev_sev*stable_stage[5]),NA))))
  
  # Matrix to store population structure through time with a row for each age and column for each time step
  pop_structure = as.data.frame(matrix(nrow=nstages, ncol=time))
  rownames(pop_structure) = 1:5
  colnames(pop_structure) = 1:ncol(pop_structure)
  
  # Add the initial population structure into the first column 
  pop_structure[,1] = P0
  
  # Create a vector to fill to be filled in with total population size
  total_pop = rep(0, times=time)
  
  for (i in 2:time) {
    pop_structure[,i] = round(whitePine_matrix %*% pop_structure[,i-1],0)
    total_pop[i]= round(sum(pop_structure[,i]),0)
  }
    
  # biomass eq from Jenkins et al. 
  biomass = as.data.frame(matrix(nrow=nstages, ncol=time))
  
  avg_dbh = c(1,5,15,30,50)
  
  for(i in 1:length(avg_dbh)){
    for(j in 1:time){
      biomass[,j] = pop_structure[,j]*exp(-2.5356 + 2.2595*log(avg_dbh[i]))
    }
  }
  
  total_carbon = rep(0, times=time)
  for (i in 1:time) {
    total_carbon[i]=sum(biomass[,i])*carbon_coeff
  }
  
  lambda = popbio::lambda(whitePine_matrix)
  
  return(list(pop_structure,total_pop,biomass,total_carbon,lambda))
}
