





forest_carbon_function = function(climate, impact, years, time_step, 
                                  initial_carbon, clim_coef,
                                  btl_r0, btl_pop0, btl_K0, btl_parms,
                                  pine_growth, pine_surv, pine_fert,
                                  pine_pop_initial = c(107,77,300,14,1)){
  
  
  clim_df = clim_df_function(years, climate)
  
  carbon_df = data.frame(year = 2010:2100, 
                        btl_r_coef = ifelse((log(sqrt(abs(clim_df$W_tmin[1]-clim_df$W_tmin)),10)+1)<1,1,
                                            (log(sqrt(abs(clim_df$W_tmin[1]-clim_df$W_tmin)),10)+1)),
                        # The btl_r_coef is used to change the r_max value for beetles each year. 
                        # Ideally this is a function of increasing minimum winter temperature because larvae mortality is also decreasing. 
                        btl_r = rep(NA,years), # Will be filled in with initial_r_max * btl_r_coef. Used to caculate the beetle population size in each time step
                        btl_pop = rep(NA,years), # Beetle population
                        p_p0 = c(pine_surv[1]*clim_coef), # Climate impact
                        p_p1 = rep(pine_surv[2], years), 
                        p_p2 = rep(pine_surv[3], years),
                        p_p3 = rep(pine_surv[4], years), # No Beetle impact (in for loop). When impacted by beetles, this value decreases as beetles increase.
                        p_p4 = rep(pine_surv[5], years), # No Beetle impact (in for loop) 
                        g01 = c(pine_growth[1]*clim_coef), # Climate impact. This value decreases as CWD increases. More water stress = less growth. 
                        g12 = c(pine_growth[2]*clim_coef), # Climate impact
                        g23 = c(pine_growth[3]*clim_coef), # Climate impact
                        g34 = c(pine_growth[4]*clim_coef), # Climate impact
                        pine0 = rep(NA,years), # Number of pines in stage 0
                        pine1 = rep(NA,years), # Number of pines in stage 1
                        pine2 = rep(NA,years), # Number of pines in stage 2
                        pine3 = rep(NA,years), # Number of pines in stage 3
                        pine4 = rep(NA,years), # Number of pines in stage 4
                        total_pine = rep(0,years), # Total pine population
                        forest_carbon = rep(0,years), # Total forest carbon, calculated more above ground biomass. 
                        fire_prob = rep(NA,years), # Fire probability, binary (High/Low)
                        fire_sev = rep(NA,years), # Fire severity, (Low, Moderate, Severe)
                        impact = rep(impact, years),
                        climate = rep(climate, years))
  
  carbon_df[,3] = round(btl_r0*carbon_df$btl_r_coef,2) # Filling in btl_r row
  carbon_df[1,4] = btl_pop0 # Initial beetle pop size
  carbon_df[1,5:9] = pine_surv # Initial pine survival
  carbon_df[1,14:18] = pine_pop_initial # Initial pine population structure
  carbon_df[1,19] = sum(carbon_df[1,14:18]) # Initial total pine population
  carbon_df[1,20] = initial_carbon # Initial carbon contained in the whitepine pop.
  carbon_df[1,21] = ifelse(clim_df$Su_tmax[1] >36 && clim_df$W_tmin[1] > 2, "H", "L") # Initial probability of fire
  carbon_df[1,22] = ifelse(clim_df$apr_snow[1] >=80, "Low", ifelse(clim_df$apr_snow[1] >=60, "Mod", "Sev")) # Initial severity of fire
  
  for(i in 2:years){
    
    ################################################# Fire ######################
    carbon_df$fire_prob[i] = ifelse(clim_df$Su_tmax[i] >36 && clim_df$W_tmin[i] > 2, "H", "L")
    carbon_df$fire_sev[i] =  ifelse(clim_df$apr_snow[i] >=80, "Low", ifelse(clim_df$apr_snow[i] >=60, "Mod", "Sev"))
    
    ################################################# Beetles ######################
    
    tmp = data.frame(time = 1:time) # Temp. data.frame for beetle output. 
    Ki = K_0*carbon_df$total_pine[i-1]/carbon_df$total_pine[1]
    beetle_parms = list(p0=btl_pop_0, r=r_0, K=Ki)
    tmp_result = ode(carbon_df$btl_pop[i-1], tmp$time, beetle_pop, beetle_parms) # ODE solver to determine beetle pop size in time 2
    carbon_df$btl_pop[i] = tmp_result[2,2] # Put ODE results into main 'carbon_df' table
    
    ################################################# Pines ######################
    
    pine_surv = as.numeric(carbon_df[i-1,5:9]) # Pine survival vector in a given time step
    pine_growth = as.numeric(carbon_df[i-1,10:13]) # Pine growth vector for a given time step
    init_pop = as.numeric(carbon_df[i-1,14:18]) # Initial population size, relative to each time step. 
    carbon_df$p_p3[i] = p3 # No beetle impact
    carbon_df$p_p4[i] = p4 # No beetle impact
    
    fire_prob = carbon_df[i,21] # Fire probability in current time step
    fire_sev = carbon_df[i,22] # Fire severity in current time step
    
    # Implement whitepine matrix model. C
    pine_fire = whitePine_matrix_FIREmodel(fertility = pine_fert,
                                           survival = pine_surv,
                                           growth = pine_growth,
                                           fire_prob = fire_prob,
                                           fire_sev = fire_sev,
                                           time= time_step,
                                           initial_pop = ntrees_0) 
    
    pine_pop = pine_fire[[1]] # Pull out pop structure
    carbon_df[i,14:18] = pine_pop[,2] # Put into main output dataframe
    carbon_df[i,19] = sum(carbon_df[i,14:18]) # Sum for the total pop size
    
    frst_crbn = pine_fire[[4]] # Pull out forest carbon 
    carbon_df[i,20] = frst_crbn[2] # Put into main output data.frame. 
    
  }
  
results(list(carbon_df))
  
}