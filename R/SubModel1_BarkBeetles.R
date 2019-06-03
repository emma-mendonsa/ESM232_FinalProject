#' Bark Beetle Logistic Growth Model
#' Add in params and other info
#' 
#' 
#' @param time time since start
#' @param P population
#' @param parms - as list with two values, r, K
#' @param r intrinsic growth rate
#' @param K carrying capacity
#' @return derivative of population with time

beetle_pop = function(Time, B, parms){

  dBeetles = parms$r * B * (1-B/parms$K)
  
  return(list(dBeetles))

  }

