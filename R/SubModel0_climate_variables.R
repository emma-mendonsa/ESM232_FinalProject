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

climate_variables_fun = function(raster,crop_layer,parms){
  
  # Crop raster using crop layer
  clim_crop = crop(raster,crop_layer)
  
  # Clim cell stats
  clim_avg = cellStats(clim_crop,"mean")
  
  # Time increments to divide channge in climate variable over
  increment = parms$total_change/parms$time
  
  # Create a slightly randomized increment vector using the above increment as the mean
  increment_vector = rnorm(mean=increment,sd=abs(increment),n=parms$time)
  
  # Create empty climate vector to store results in
  clim_vector = rep(NA,parms$time)
  
  # Initial value of climate variables
  clim_vector[1] = clim_avg
  
  # Fill in climate vector
  for(i in 2:parms$time){
    clim_vector[i] = clim_vector[i-1]+increment_vector[i]  
  }
  
  return(clim_vector)
  
}
