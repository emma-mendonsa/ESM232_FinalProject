#' Whitebark pine projection model
#'
#' @param stages names of each stage to be used to name rows and columns in matrix
#' @param fertility vector of fertility values for each stage
#' @param survival vector of survival values for each stage. Values represent the probability of surviving and remaining in the same stage
#' @param growth vector of growth values for each stage. Values represent the probability of surviving and growing into the next stage
#' @param initial_pop vector of the initial population struture
#' @param time length of time to calculate the population structure over
#' @return List that contains the population structure through time, total population through time, and the population lambda (finite growth rate) based on the initial matrix
#' @references Jules et al. (2016) The relative contributions of disease and insects in the decline of a long-lived tree: stochastic demographic model of whitebark pine (Pinus albicaulis). Forest Ecology and Managment. 381. Pages 144-156. 

whitePine_matrix_model = function(stages, fertility, survival, growth, 
                                  initial_pop, time) {
  
  nstages = length(fertility)
  
  ### Error checking to ensure the right data exists for each stage###
  # Survivability value for each stae
  if ((nstages!=length(survival) ))
  {return("Missing or extra data: Length of fertility doesn’t match survivability")}
  
  # Initialize Leslie Matrix
  whitePine_matrix = matrix(nrow=nstages, ncol=nstages)
  colnames(whitePine_matrix) = stages
  rownames(whitePine_matrix) = stages
  
  # Fill in with zero values
  whitePine_matrix[,] = 0.0
  
  # Add fertility values to matric
  whitePine_matrix[1,] = fertility
  
  # Add survival values to matrix 
  for (i in 1:(nstages)) {
    whitePine_matrix[i,i] = survival[i]
  }
 
  # Add growth values to matrix
  for (i in 1:(nstages-1)) {
    whitePine_matrix[i+1,i] = growth[i]
  } 
  
  # Matrix to store population structure through time with a row for each age and column for each time step
  pop_structure = as.data.frame(matrix(nrow=nstages, ncol=time))
  rownames(pop_structure) = 1:5
  colnames(pop_structure) = 1:ncol(pop_structure)
  
  # Add the initial population structure into the first column 
  pop_structure[,1] = initial_pop
  
  # Create a vector to be filled in with total population size
  total_pop = rep(0, times=time)
  
  for (i in 2:time) {
    # Calculate whitebark pine population structure through time
    pop_structure[,i] = round(whitePine_matrix %*% pop_structure[,i-1],0)
    
    # Find total population size in each timestep
    total_pop[i]=round(sum(pop_structure[,i]),0)
  }
  
  # Find the population growth rate of the initial matrix
  lambda = popbio::lambda(whitePine_matrix)
  
  return(list(pop_structure=pop_structure, # Table of the population structure
              total_pop=total_pop, # Total population size
              lambda=lambda)) # Population growth rate
}
