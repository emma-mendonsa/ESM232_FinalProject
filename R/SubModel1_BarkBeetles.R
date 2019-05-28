#' Bark Beetle Logistic Growth Model
#' 
#' @param r_beetle Beetle growth rate

beetle_pop = function(N,r,K){
  
  N = N + r*N*(1-N/K)
  
  return(N)
  
}
