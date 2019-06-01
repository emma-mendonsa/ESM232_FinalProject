climate_change_fun = function(total_change,time,init_value){
  increment = total_change/time
  increment_vector = rnorm(mean=increment,sd=abs(increment),n=time)
  
  clim_vector = rep(NA,time)
  clim_vector[1] = init_value
  
  for(i in 2:time){
    clim_vector[i] = clim_vector[i-1]+increment_vector[i]  
  }
  
  return(clim_vector)
  
}
