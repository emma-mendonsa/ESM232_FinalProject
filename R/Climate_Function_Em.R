

clim_df_function = function(years, climate){
  
  source("../R/SubModel0_climate_variables.R")
  
  
  # Find average tmin and tmax values for Sierra Nevada Ecoregion based on historic data
  cwd_raw = raster("../climate_data/cwd1981_2010_ave_HST_1559416632/cwd1981_2010_ave_HST_1559416632.tif")
  snwpck_raw = raster("../climate_data/aprpck1981_2010_ave_HST_1559416720/aprpck1981_2010_ave_HST_1559416720.tif")
  W_tmin_raw = raster("../climate_data/tmn1981_2010djf_ave_HST_1559415579/tmn1981_2010djf_ave_HST_1559415579.tif")
  Sp_tmin_raw = raster("../climate_data/tmn1981_2010mar_ave_HST_1559415612/tmn1981_2010mar_ave_HST_1559415612.tif")
  F_tmin_raw = raster("../climate_data/tmn1981_2010nov_ave_HST_1559415731/tmn1981_2010nov_ave_HST_1559415731.tif")
  Su_tmax_raw = raster("../climate_data/tmx1981_2010jja_ave_HST_1559415851/tmx1981_2010jja_ave_HST_1559415851.tif")
  crop_layer = readOGR("../climate_data/SNC_Boundary_Shapefile/SNC_Boundary.shp")
  

  #Create two climate scenarios:

  ##Business as Usual Climate Dataframe##
  clim_BAU = data.frame(year = 2010:2100,
                        time = 1:years)
  
  # Fill in climate data table
  clim_parms = list(total_change=-50,time=years)
  clim_BAU$apr_snow = climate_variables_fun(snwpck_raw,crop_layer = crop_layer,parms = clim_parms)
  
  clim_parms = list(total_change=400,time=years)
  clim_BAU$cwd = climate_variables_fun(cwd_raw,crop_layer = crop_layer,parms = clim_parms)
  
  clim_parms = list(total_change=10,time=years)
  clim_BAU$F_tmin = climate_variables_fun(F_tmin_raw,crop_layer = crop_layer,parms = clim_parms)
  
  clim_parms = list(total_change=10,time=years)
  clim_BAU$Sp_tmin = climate_variables_fun(Sp_tmin_raw,crop_layer = crop_layer,parms = clim_parms)
  
  clim_parms = list(total_change=10,time=years)
  clim_BAU$W_tmin = climate_variables_fun(W_tmin_raw,crop_layer = crop_layer,parms = clim_parms)
  
  clim_parms = list(total_change=10,time=years)
  clim_BAU$Su_tmax = climate_variables_fun(Su_tmax_raw,crop_layer = crop_layer,parms = clim_parms)
  
 
  
  ##Moderate Climate Dataframe##
  clim_moderate = data.frame(year = 2010:2100,
                       time = 1:years,
                       apr_snow = rep(NA,years),
                       cwd = rep(NA,years),
                       F_tmin = rep(NA,years),
                       W_tmin = rep(NA,years),
                       Sp_tmin = rep(NA,years),
                       Su_tmax = rep(NA,years))
  
  
  # Fill in climate data table
  clim_parms = list(total_change=-20,time=years)
  clim_moderate$apr_snow = climate_variables_fun(snwpck_raw,crop_layer = crop_layer,parms = clim_parms)
  
  clim_parms = list(total_change=200,time=years)
  clim_moderate$cwd = climate_variables_fun(cwd_raw,crop_layer = crop_layer,parms = clim_parms)
  
  clim_parms = list(total_change=5,time=years)
  clim_moderate$F_tmin = climate_variables_fun(F_tmin_raw,crop_layer = crop_layer,parms = clim_parms)
  
  clim_parms = list(total_change=5,time=years)
  clim_moderate$Sp_tmin = climate_variables_fun(Sp_tmin_raw,crop_layer = crop_layer,parms = clim_parms)
  
  clim_parms = list(total_change=5,time=years)
  clim_moderate$W_tmin = climate_variables_fun(W_tmin_raw,crop_layer = crop_layer,parms = clim_parms)
  
  clim_parms = list(total_change=5,time=years)
  clim_moderate$Su_tmax = climate_variables_fun(Su_tmax_raw,crop_layer = crop_layer,parms = clim_parms)




  
  clim_df = ifelse(climate == "BAU", clim_BAU, 
            ifelse(climate == "moderate", clim_moderate,
                   "Error: Invalid climate entry"))
  
  return(clim_df)
  
}

