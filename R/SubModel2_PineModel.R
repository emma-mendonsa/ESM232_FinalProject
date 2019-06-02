#'
#'

whitePine_matrix_model = function(fertility, survival, growth, initial_pop, time, carbon_coeff = 0.47) {
  
  nstages = length(fertility)
  
  ### Error checking to ensure the right data exists for each stage###
  # Survivability value for each stae
  if ((nstages!=length(survival) ))
  {return("Missing or extra data: Length of fertility doesn’t match survivability")}
  
  if ((nstages!=length(initial_pop) ))
  { return("Missing or extra data: Length of initial population doesn’t match length of fertility") }
  
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
  
  lambda = popbio::lambda(whitePine_matrix)
  
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
    total_pop[i]=round(sum(pop_structure[,i]),0)
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
  
  return(list(pop_structure,total_pop,biomass,total_carbon))
}
