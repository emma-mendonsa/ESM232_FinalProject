#' Whitebark pine population dynamics with influences from forest fires, climate change, and mountain pine beetles. 
#' 
#' This function brings together submodels 1 through 4 which calculate the mountain pine beetle population size, stage structure matrix population model for whitebark pines, forest fire probability and severity, and the above ground biomass and carbon in the forest of interest, respectively. Notes throughout will explain how the submodels interact. 
#' 
#' @param clim_df Dataframe of climate variables through time, probuced using submodel 0 
#' @param clim_scen_name Name of the climate change scenario. This is used in the output data.frame to specify the climate impact 
#' @param impact_type Dynamics having direct impacts on the whitebark pine population and beetle population.
#' @param scenario_parms List containing the total length of time being analyzed (i.e. 91 years) and the timestep of the analysis (i.e. two years)
#' @param beetle_parms List of all input variables for the mountain pine beetle logistic growth model. See submodel 1 for details
#' @param pine_parms List of all input variables for the whitebark pine matrix projection model. See submodel 2 for details
#' @param beetle_K beetle carrying capacity in the initial time step/ 
#' 
#' @return output dataframe. See where output_df is initialized below for specific attributes. 
#' @author Emma Mendonsa and Claire Powers
#'  

whitePine_PopDynamics_model = function(clim_df,clim_scen_name,impact_type,scenario_parms,pine_parms,beetle_parms,beetle_K){
  
  ########## Climate only impact ############
  if(impact_type=="climate_only"){
  
  clim_coef = clim_df$cwd[1]/clim_df$cwd # Coefficient used to proportionally change different whitebark pine vital rates as the climate changes. Increases in CWD decrease seedling survival (the youngest age class) and decrease the growth rates. 
    
  years = nrow(clim_df) # Determine number of years based on number of rows in the climate data.frame. 
  
  output_df = data.frame(year = clim_df$year, # Year column
                         # Increase the beetle intrinsic growth rate as minimum winter temperatures increase
                        btl_rmax_coef = ifelse((log(sqrt(abs(clim_df$W_tmin[1]-clim_df$W_tmin)),10)+1)<1,1,
                                            (log(sqrt(abs(clim_df$W_tmin[1]-clim_df$W_tmin)),10)+1)),
                        btl_rmax = rep(NA,years), # Will be filled in with initial_r_max * btl_r_coef. 
                        btl_pop = rep(NA,years), # Beetle population
                        p_p0 = c (pine_parms$pine_surv[1]*clim_coef), # Climate impact. More water stress = fewer seedlings survive
                        p_p1 = rep(pine_parms$pine_surv[2], years), 
                        p_p2 = rep(pine_parms$pine_surv[3], years),
                        p_p3 = rep(pine_parms$pine_surv[4], years), # No Beetle impact. When impacted by beetles, this value decreases as beetles increase.
                        p_p4 = rep(pine_parms$pine_surv[5], years), # No Beetle impact.
                        g01 = c(pine_parms$pine_growth[1]*clim_coef), # Climate impact. This value decreases as CWD increases. More water stress = less growth. 
                        g12 = c(pine_parms$pine_growth[1]*clim_coef), # Climate impact
                        g23 = c(pine_parms$pine_growth[1]*clim_coef), # Climate impact
                        g34 = c(pine_parms$pine_growth[1]*clim_coef), # Climate impact
                        pine0 = rep(NA,years), # Number of pines in stage 0
                        pine1 = rep(NA,years), # Number of pines in stage 1
                        pine2 = rep(NA,years), # Number of pines in stage 2
                        pine3 = rep(NA,years), # Number of pines in stage 3
                        pine4 = rep(NA,years), # Number of pines in stage 4
                        total_pine = rep(0,years), # Total pine population
                        preFire_forest_carbon = rep(0,years), # Total forest carbon (kg), calculated from above ground biomass(kg). 
                        fire_prob = rep(NA,years), # Fire probability, binary (High/Low)
                        fire_sev = rep(NA,years), # Fire severity, (Low, Moderate, Severe)
                        postFire_forest_carbon = rep(0,years), # Total forest carbon (kg) after a fire
                        impact_type = rep(impact_type,years), # Impact type
                        climate = rep(clim_scen_name, years)) # Climate change scenario
  
  # This section places intial values into the data.frame created above. 
  output_df[,3] = round(beetle_parms$r_max*output_df$btl_rmax_coef,2) # Beetle intrinsic growth rate through time. Increases as minimum winter temperatures increase.
  output_df[1,4] = beetle_parms$btl_pop_0 # Initial beetle pop size
  output_df[1,5:9] = pine_parms$pine_surv # Initial pine survival rates
  output_df[1,14:18] = pine_parms$initial_pop # Initial pine population structure
  output_df[1,19] = sum(output_df[1,14:18]) # Initial total pine population
  output_df[1,20] = pine_parms$initial_carbon # Initial carbon contained in the whitepine population.
  output_df[1,21] = ifelse(clim_df$Su_tmax[1] >36 && clim_df$W_tmin[1] > 2, "H", "L") # Initial probability of fire
  output_df[1,22] = ifelse(clim_df$apr_snow[1] >=80, "Low", ifelse(clim_df$apr_snow[1] >=60, "Mod", "Sev")) # Initial severity of fire
  output_df[1,23] = pine_parms$initial_carbon # Initial carbon contained in the whitepine population.
  
  for(i in 2:nrow(clim_df)){
    
    ################################# Beetles ######################################
    # Calculate beetle population dynamics using submodel 1
    
    tmp = data.frame(time = 1:scenario_parms$time_step) # Temporary data.frame to be used by the ode() function
    
    # Update beetle carrying capacity based on size of pine population in previous and initial time step. 
    # This causes carrying capacity to change proportionally with the size of the forest
    Ki = ifelse(output_df$total_pine[i-1]>0, # If the total pine population is greater than 0
                beetle_K*output_df$total_pine[i-1]/output_df$total_pine[1], # The the beetle carrying capacity increases
                beetle_K) #Otherwise it is the initial carrying capacity
    
    # Reset beetle parameters with new rmax and carrying capacity
    beetle_parms = list(btl_pop_0=beetle_parms$btl_pop_0, 
                        r_max = output_df$btl_rmax[i], 
                        K_0_btl=Ki)
    
    tmp_result = ode(beetle_parms$btl_pop_0, tmp$time, beetle_pop, beetle_parms) # ODE solver to determine beetle pop size in time 2
    output_df$btl_pop[i] = round(tmp_result[2,2],0)  # Put ODE results into main 'output_df' table
    
    ################################################### Pines #####################################################
    # Calculate pine population dynamics using submodel 2
    
    pine_surv = as.numeric(output_df[i-1,5:9]) # Pine survival vector in a given time step
    pine_growth = as.numeric(output_df[i-1,10:13]) # Pine growth vector for a given time step
    init_pop = as.numeric(output_df[i-1,14:18]) # Initial population size, relative to each time step. 
    output_df$p_p3[i] = p3 # No beetle impact
    output_df$p_p4[i] = p4 # No beetle impact
    
    # Implement whitebark pine matrix model
    pine = whitePine_matrix_model(stages = pine_parms$stages,
                                  fertility = pine_parms$pine_fert,
                                  survival = pine_surv,
                                  growth = pine_growth,
                                  time= scenario_parms$time_step,
                                  initial_pop = init_pop) 
    
    output_df[i,14:18] = pine$pop_structure[,2] # Put pop structure into main output dataframe
    output_df[i,19] = pine$total_pop[2] # Put total pop into main output dataframe
    
    ################################################ Pre-fire forest carbon ###########################################
    # Calculate forest carbon using submodel 4
    
    frst_crbn = calc_biomass_carbon(avg_stage_dbh = pine_parms$avg_stage_dbh, pop_structure = as.numeric(output_df[i,14:18]))
    output_df[i,20] = frst_crbn[2] # Put into main output data.frame. 
    
    ################################################# Fire ##################################################
    # Calculate the fire dynamics and impacts using submodel 3
    
    output_df$fire_prob[i] = ifelse(clim_df$Su_tmax[i] > 36 && clim_df$W_tmin[i] > 2, "H", "L") # Determine fire probabilty based on summer and winter temperatures
    output_df$fire_sev[i] = ifelse(clim_df$apr_snow[i] >= 80, "Low", ifelse(clim_df$apr_snow[i] >=60, "Mod", "Sev")) # Determine fire severity based on april snowpack
    
    fire = fire_model(initial_pop = as.numeric(output_df[i,14:18]),
                      pine_biomass = frst_crbn$total_biomass,
                      fire_prob = output_df$fire_prob[i],
                      fire_sev = output_df$fire_sev[i],
                      forest_carbon = frst_crbn$total_carbon)
    
    output_df[i,14:18] = fire$postFire_pop
    output_df$postFire_forest_carbon[i] = fire$post_fire_carbon
  
    } # End for loop
  return(output_df) 
  
  }else  
    if(impact_type=="beetles_only"){ ########## Beetle only impact - Notes are placed where values/dynamics differ from the climate only impact
      
      clim_coef = 1 # No climate impact so CWD does not effect growth or survival rates
      years = nrow(clim_df)
      
      output_df = data.frame(year = clim_df$year, 
                             btl_rmax_coef = ifelse((log(sqrt(abs(clim_df$W_tmin[1]-clim_df$W_tmin)),10)+1)<1,1,
                                                    (log(sqrt(abs(clim_df$W_tmin[1]-clim_df$W_tmin)),10)+1)), 
                             btl_rmax = rep(NA,years), 
                             btl_pop = rep(NA,years), 
                             p_p0 = c (pine_parms$pine_surv[1]* clim_coef), # Clim_coef = 1 so no impact
                             p_p1 = rep(pine_parms$pine_surv[2],years), 
                             p_p2 = rep(pine_parms$pine_surv[3],years),
                             p_p3 = rep(pine_parms$pine_surv[4],years), 
                             p_p4 = rep(pine_parms$pine_surv[5],years), 
                             g01 = c(pine_parms$pine_growth[1]*clim_coef), # Clim_coef = 1 so no impact 
                             g12 = c(pine_parms$pine_growth[1]*clim_coef), # Clim_coef = 1 so no impact 
                             g23 = c(pine_parms$pine_growth[1]*clim_coef), # Clim_coef = 1 so no impact 
                             g34 = c(pine_parms$pine_growth[1]*clim_coef), # Clim_coef = 1 so no impact 
                             pine0 = rep(NA,years),
                             pine1 = rep(NA,years),
                             pine2 = rep(NA,years), 
                             pine3 = rep(NA,years), 
                             pine4 = rep(NA,years), 
                             total_pine = rep(0,years), 
                             preFire_forest_carbon = rep(0,years), 
                             fire_prob = rep(NA,years), 
                             fire_sev = rep(NA,years), 
                             postFire_forest_carbon = rep(0,years),
                             impact_type = rep(impact_type,years),
                             climate = rep(clim_scen_name, years))
      
      # This section places intial values into the data.frame created above. 
      output_df[,3] = round(beetle_parms$r_max*output_df$btl_rmax_coef,2)
      output_df[1,4] = beetle_parms$btl_pop_0 
      output_df[1,5:9] = pine_parms$pine_surv 
      output_df[1,14:18] = pine_parms$initial_pop 
      output_df[1,19] = sum(output_df[1,14:18]) 
      output_df[1,20] = pine_parms$initial_carbon 
      output_df[1,21] = ifelse(clim_df$Su_tmax[1] >36 && clim_df$W_tmin[1] > 2, "H", "L") 
      output_df[1,22] = ifelse(clim_df$apr_snow[1] >=80, "Low", ifelse(clim_df$apr_snow[1] >=60, "Mod", "Sev")) 
      output_df[1,23] = pine_parms$initial_carbon
      
      for(i in 2:91){
        
        ################################# Beetles ######################################
        
        tmp = data.frame(time = 1:scenario_parms$time_step) 
        Ki = ifelse(output_df$total_pine[i-1]>0,beetle_K*output_df$total_pine[i-1]/output_df$total_pine[1],beetle_K)
        
        beetle_parms = list(btl_pop_0=beetle_parms$btl_pop_0, 
                            r_max = output_df$btl_rmax[i], 
                            K_0_btl=Ki)
        
        tmp_result = ode(beetle_parms$btl_pop_0, tmp$time, beetle_pop, beetle_parms) 
        output_df$btl_pop[i] = round(tmp_result[2,2],0)  
        
        ################################### Pines #######################################
        
        pine_surv = as.numeric(output_df[i-1,5:9]) 
        pine_growth = as.numeric(output_df[i-1,10:13]) 
        init_pop = as.numeric(output_df[i-1,14:18]) 
        
        # As the beetle population increases, the survival of stage 3 and 4 whitebark pines decrease
        output_df$p_p3[i] = ifelse(output_df$btl_pop[i]>1,output_df$p_p3[i-1]*log(output_df$btl_pop[1])/log(output_df$btl_pop[i]),output_df$p_p3[1]) # beetle impact
        output_df$p_p4[i] = ifelse(output_df$btl_pop[i]>1,output_df$p_p4[i-1]*log(output_df$btl_pop[1])/log(output_df$btl_pop[i]),output_df$p_p4[1]) # beetle impact
        
        
        # Implement whitebark pine matrix model
        pine = whitePine_matrix_model(stages = pine_parms$stages,
                                      fertility = pine_parms$pine_fert,
                                      survival = pine_surv,
                                      growth = pine_growth,
                                      time= scenario_parms$time_step,
                                      initial_pop = init_pop) 
        
        output_df[i,14:18] = pine$pop_structure[,2] 
        output_df[i,19] = pine$total_pop[2] 
        
        ################################################ Pre-fire forest carbon ################################################
        
        frst_crbn = calc_biomass_carbon(avg_stage_dbh = pine_parms$avg_stage_dbh, pop_structure = as.numeric(output_df[i,14:18]))
        output_df[i,20] = frst_crbn[2] 
        
        ################################################# Fire ##################################################
        
        output_df$fire_prob[i] = ifelse(clim_df$Su_tmax[i] > 36 && clim_df$W_tmin[i] > 2, "H", "L")
        output_df$fire_sev[i] = ifelse(clim_df$apr_snow[i] >= 80, "Low", ifelse(clim_df$apr_snow[i] >=60, "Mod", "Sev"))
        
        fire = fire_model(initial_pop = as.numeric(output_df[i,14:18]),
                          pine_biomass = frst_crbn$total_biomass,
                          fire_prob = output_df$fire_prob[i],
                          fire_sev = output_df$fire_sev[i],
                          forest_carbon = frst_crbn$total_carbon)
        
        output_df[i,14:18] = fire$postFire_pop
        output_df$postFire_forest_carbon[i] = fire$post_fire_carbon
      } # End for loop
      
      return(output_df) 
       
    }
  
  else
    if(impact_type == "beetles_and_climate"){ ###### Beetles and climate impact ######### Combination of beetle only and climate only impacts
      
      clim_coef = clim_df$cwd[1]/clim_df$cwd
      
      years = nrow(clim_df)
      output_df = data.frame(year = clim_df$year, 
                             btl_rmax_coef = ifelse((log(sqrt(abs(clim_df$W_tmin[1]-clim_df$W_tmin)),10)+1)<1,1,
                                                    (log(sqrt(abs(clim_df$W_tmin[1]-clim_df$W_tmin)),10)+1)),
                             btl_rmax = rep(NA,years),
                             btl_pop = rep(NA,years),
                             p_p0 = c (pine_parms$pine_surv[1]*  clim_coef), # Climate impact
                             p_p1 = rep(pine_parms$pine_surv[2], years), 
                             p_p2 = rep(pine_parms$pine_surv[3], years),
                             p_p3 = rep(pine_parms$pine_surv[4], years), 
                             p_p4 = rep(pine_parms$pine_surv[5], years),
                             g01 = c(pine_parms$pine_growth[1]* clim_coef), # Climate impact
                             g12 = c(pine_parms$pine_growth[1]* clim_coef), # Climate impact
                             g23 = c(pine_parms$pine_growth[1]* clim_coef), # Climate impact
                             g34 = c(pine_parms$pine_growth[1]* clim_coef), # Climate impact
                             pine0 = rep(NA,years), 
                             pine1 = rep(NA,years), 
                             pine2 = rep(NA,years),
                             pine3 = rep(NA,years), 
                             pine4 = rep(NA,years), 
                             total_pine = rep(0,years),
                             preFire_forest_carbon = rep(0,years), 
                             fire_prob = rep(NA,years),
                             fire_sev = rep(NA,years),
                             postFire_forest_carbon = rep(0,years),
                             impact_type = rep(impact_type,years),
                             climate = rep(clim_scen_name, years))
      
      # This section places intial values into the data.frame created above. 
      output_df[,3] = round(beetle_parms$r_max*output_df$btl_rmax_coef,2) 
      output_df[1,4] = beetle_parms$btl_pop_0 
      output_df[1,5:9] = pine_parms$pine_surv 
      output_df[1,14:18] = pine_parms$initial_pop 
      output_df[1,19] = sum(output_df[1,14:18]) 
      output_df[1,20] = pine_parms$initial_carbon 
      output_df[1,21] = ifelse(clim_df$Su_tmax[1] >36 && clim_df$W_tmin[1] > 2, "H", "L") 
      output_df[1,22] = ifelse(clim_df$apr_snow[1] >=80, "Low", ifelse(clim_df$apr_snow[1] >=60, "Mod", "Sev")) 
      output_df[1,23] = pine_parms$initial_carbon
      
      for(i in 2:nrow(clim_df)){
        
        ################################# Beetles ######################################
        
        tmp = data.frame(time = 1:scenario_parms$time_step) 
        
        Ki = ifelse(output_df$total_pine[i-1]>0,beetle_K*output_df$total_pine[i-1]/output_df$total_pine[1],beetle_K)
        
        # Reset beetle parameters with new rmax and carrying capacity
        beetle_parms = list(btl_pop_0=beetle_parms$btl_pop_0, 
                            r_max = output_df$btl_rmax[i], 
                            K_0_btl=Ki)
        
        tmp_result = ode(beetle_parms$btl_pop_0, tmp$time, beetle_pop, beetle_parms)
        output_df$btl_pop[i] = round(tmp_result[2,2],0) 
        
        ################################### Pines #######################################
        
        pine_surv = as.numeric(output_df[i-1,5:9]) 
        pine_growth = as.numeric(output_df[i-1,10:13]) 
        init_pop = as.numeric(output_df[i-1,14:18]) 
        output_df$p_p3[i] = ifelse(output_df$btl_pop[i]>1,output_df$p_p3[i-1]*log(output_df$btl_pop[1])/log(output_df$btl_pop[i]),output_df$p_p3[1]) # beetle impact
        output_df$p_p4[i] = ifelse(output_df$btl_pop[i]>1,output_df$p_p4[i-1]*log(output_df$btl_pop[1])/log(output_df$btl_pop[i]),output_df$p_p4[1]) # beetle impact
        
        # Implement whitebark pine matrix model
        pine = whitePine_matrix_model(stages = pine_parms$stages,
                                      fertility = pine_parms$pine_fert,
                                      survival = pine_surv,
                                      growth = pine_growth,
                                      time= scenario_parms$time_step,
                                      initial_pop = init_pop) 
        
        output_df[i,14:18] = pine$pop_structure[,2] 
        output_df[i,19] = pine$total_pop[2] 
        
        ################################################ Pre-fire forest carbon ################################################
        
        frst_crbn = calc_biomass_carbon(avg_stage_dbh = pine_parms$avg_stage_dbh, pop_structure = as.numeric(output_df[i,14:18]))
        output_df[i,20] = frst_crbn[2] 
        
        ################################################# Fire ##################################################
        
        output_df$fire_prob[i] = ifelse(clim_df$Su_tmax[i] > 36 && clim_df$W_tmin[i] > 2, "H", "L")
        output_df$fire_sev[i] = ifelse(clim_df$apr_snow[i] >= 80, "Low", ifelse(clim_df$apr_snow[i] >=60, "Mod", "Sev"))
        
        fire = fire_model(initial_pop = as.numeric(output_df[i,14:18]),
                          pine_biomass = frst_crbn$total_biomass,
                          fire_prob = output_df$fire_prob[i],
                          fire_sev = output_df$fire_sev[i],
                          forest_carbon = frst_crbn$total_carbon)
        
        output_df[i,14:18] = fire$postFire_pop
        output_df$postFire_forest_carbon[i] = fire$post_fire_carbon
      } # End for loop
      return(output_df)
      
    }
}


