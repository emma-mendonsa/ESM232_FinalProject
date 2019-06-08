#' Wrapper function
#' 
#' 
#' 

whitePine_PopDynamics_model = function(clim_df,clim_scen_name,scenario_parms,pine_parms,beetle_parms){
  
  years = nrow(clim_df)
  
  output_df = data.frame(year = clim_df$year, 
                        btl_rmax_coef = ifelse((log(sqrt(abs(clim_df$W_tmin[1]-clim_df$W_tmin)),10)+1)<1,1,
                                            (log(sqrt(abs(clim_df$W_tmin[1]-clim_df$W_tmin)),10)+1)),
                        btl_r = rep(NA,years), # Will be filled in with initial_r_max * btl_r_coef. Used to caculate the beetle population size in each time step
                        btl_pop = rep(NA,years), # Beetle population
                        p_p0 = c (pine_parms$pine_surv[1]* scenario_parms$clim_coef), # Climate impact
                        p_p1 = rep(pine_parms$pine_surv[2], years), 
                        p_p2 = rep(pine_parms$pine_surv[3], years),
                        p_p3 = rep(pine_parms$pine_surv[4], years), # No Beetle impact (in for loop). When impacted by beetles, this value decreases as beetles increase.
                        p_p4 = rep(pine_parms$pine_surv[5], years), # No Beetle impact (in for loop) 
                        g01 = c(pine_parms$pine_growth[1]*scenario_parms$clim_coef), # Climate impact. This value decreases as CWD increases. More water stress = less growth. 
                        g12 = c(pine_parms$pine_growth[1]*scenario_parms$clim_coef), # Climate impact
                        g23 = c(pine_parms$pine_growth[1]*scenario_parms$clim_coef), # Climate impact
                        g34 = c(pine_parms$pine_growth[1]*scenario_parms$clim_coef), # Climate impact
                        pine0 = rep(NA,years), # Number of pines in stage 0
                        pine1 = rep(NA,years), # Number of pines in stage 1
                        pine2 = rep(NA,years), # Number of pines in stage 2
                        pine3 = rep(NA,years), # Number of pines in stage 3
                        pine4 = rep(NA,years), # Number of pines in stage 4
                        total_pine = rep(0,years), # Total pine population
                        preFire_forest_carbon = rep(0,years), # Total forest carbon, calculated more above ground biomass. 
                        fire_prob = rep(NA,years), # Fire probability, binary (High/Low)
                        fire_sev = rep(NA,years), # Fire severity, (Low, Moderate, Severe)
                        postFire_forest_carbon = rep(0,years),
                        impact = rep(impact, years),
                        climate = rep(clim_scen_name, years))
  
  output_df[,3] = round(beetle_parms$r_max_0*output_df$btl_rmax_coef,2) # Filling in btl_r row
  output_df[1,4] = beetle_parms$btl_pop_0 # Initial beetle pop size
  output_df[1,5:9] = pine_parms$pine_surv # Initial pine survival
  output_df[1,14:18] = pine_parms$initial_pop # Initial pine population structure
  output_df[1,19] = sum(output_df[1,14:18]) # Initial total pine population
  output_df[1,20] = pine_parms$initial_carbon # Initial carbon contained in the whitepine pop.
  output_df[1,21] = ifelse(clim_df$Su_tmax[1] >36 && clim_df$W_tmin[1] > 2, "H", "L") # Initial probability of fire
  output_df[1,22] = ifelse(clim_df$apr_snow[1] >=80, "Low", ifelse(clim_df$apr_snow[1] >=60, "Mod", "Sev")) # Initial severity of fire
  
  for(i in 2:nrow(clim_df)){
    
    ################################# Beetles ######################################
    
    tmp = data.frame(time = 1:scenario_parms$time_step) # Temporary data.frame for beetle output. 
    Ki = beetle_parms$K_0_btl*output_df$total_pine[i-1]/output_df$total_pine[1]
    beetle_parms = list(p0=btl_pop_0, r=r_0, K=Ki)
    tmp_result = ode(output_df$btl_pop[i-1], tmp$time, beetle_pop, beetle_parms) # ODE solver to determine beetle pop size in time 2
    output_df$btl_pop[i] = tmp_result[2,2] # Put ODE results into main 'output_df' table
    
    ################################### Pines #######################################
    
    pine_surv = as.numeric(output_df[i-1,5:9]) # Pine survival vector in a given time step
    pine_growth = as.numeric(output_df[i-1,10:13]) # Pine growth vector for a given time step
    init_pop = as.numeric(output_df[i-1,14:18]) # Initial population size, relative to each time step. 
    output_df$p_p3[i] = p3 # No beetle impact
    output_df$p_p4[i] = p4 # No beetle impact
    
    # Implement whitebark pine matrix model
    pine = whitePine_matrix_model(fertility = pine_fert,
                                  survival = pine_surv,
                                  growth = pine_growth,
                                  time=2,
                                  initial_pop = ntrees_0) 
    
    output_df[i,14:18] = pine$pop_structure[,2] # Put pop structure into main output dataframe
    output_df[i,19] = pine$total_pop[2] # Put total pop into main output dataframe
    
    
    ################################################ Pre-fire forest carbon #####
    
    avg_stage_dbh = c(1,5,15,30,50)
    frst_crbn = calc_biomass_carbon(avg_stage_dbh = avg_stage_dbh, pop_structure = as.numeric(output_df[i,14:18]))
    output_df[i,20] = frst_crbn[2] # Put into main output data.frame. 
    
    ################################################# Fire ######################
    
    output_df$fire_prob[i] = ifelse(clim_BAU$Su_tmax[i] >36 && clim_BAU$W_tmin[i] > 2, "H", "L")
    output_df$fire_sev[i] =  ifelse(clim_BAU$apr_snow[i] >=80, "Low", ifelse(clim_BAU$apr_snow[i] >=60, "Mod", "Sev"))
    
    fire = burnt_biomass_model(initial_pop = as.numeric(output_df[i,14:18]),
                               pine_biomass = frst_crbn$total_biomass,
                               fire_prob = output_df$fire_prob[i],
                               fire_sev = output_df$fire_sev[i],
                               forest_carbon = frst_crbn$total_carbon)
    
    output_df[i,14:18] = fire$postFire_pop
    output_df$postFire_forest_carbon = fire$post_fire_carbon
    
  }
  
  return(output_df)
  
}