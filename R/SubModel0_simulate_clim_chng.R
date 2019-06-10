#' Simple simulation of climate change 
#' 
#' @param raster A raster layer of historic climate data
#' @param crop_layer A layer of the extent of the area of interest
#' @param parms A list of parameters that specifies the total change (i.e. 10 degrees celcius) and the time span you'd like change to occur over (i.e. 91 years)
#'
#' @return clim_vector vector of the value of the climate variable through time
#' 

simulate_clim_chng_fun = function(raster,crop_layer,parms){
  
  # Crop raster using crop layer
  clim_crop = crop(raster,crop_layer) # Crop the raster layer to the outline of the area of interest
  
  # Clim cell stats
  clim_avg = cellStats(clim_crop,"mean") # Find the mean value of that cropped raster
  
  # Time increments to divide channge in climate variable over
  increment = parms$total_change/parms$time # Find how much the mean value needs to change each time step
  
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
