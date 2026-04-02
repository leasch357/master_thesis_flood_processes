event_classification <- function(input_data, model_output,warmup, param) { #prepare decision variables for flood process classification
  #method 1 uses relative soil moisture based on fc (sat_start=sm/fc)
  #input:
  #input_data model forcing and observed flow
  #model_output hbv output containing qsim and states
  #warmup number of warmup days to remove
  #param calibrated hbv parameter set (needs fc)
  #output:
  #decision df with flood variables (indicators, thresholds)
  
  #prepare flood_data_df
  #remove warmup
  output_data <- data.frame(date = input_data$date, Q_obs = input_data$flow)
  output_data <- slice(output_data, (warmup + 1):nrow(output_data))
  #ensure series starts on jan 1 by dropping the first incomplete year
  first_year <- year(output_data$date[1])
  if (yday(output_data$date[1]) != 1) {
    output_data <- output_data %>% filter(year(date) > first_year)
  }
  
  #extract annual maximum observed discharge (amax) 
  flood_data_df <- output_data %>%
    mutate(jahr= year(date)) %>%
    group_by(jahr) %>% 
    slice_max(order_by = Q_obs, n = 1, with_ties = FALSE) %>% 
    ungroup()
  
  #keep date and magnitude and define t0 (start of 7 day window)
  flood_data_df <- flood_data_df[, c("date", "Q_obs")]
  colnames(flood_data_df) = c("flood_date", "flood_magn")
  flood_data_df$t0 <- flood_data_df$flood_date - 6 
  
  #extract hbv state variables (snow and soil moisture)
  storage_alle <- attr(model_output, "STATES")
  length_storage_df <- length(storage_alle[,1])
  n_dates  <- nrow(input_data)
  
  #build snow storages time series and remove warmup period
  if (length_storage_df == n_dates){
    snowpack <- data.frame(date=input_data$date, snow_component =   storage_alle[,1], liquid_component = storage_alle[,2])
  } else {
    snowpack <- data.frame(date=input_data$date, snow_component = storage_alle[2:length_storage_df,1], liquid_component = storage_alle[2:length_storage_df,2])
  }
  snowpack <- slice(snowpack, (warmup + 1):nrow(snowpack))
  
  #compute 7 day snowmelt (positive decreases in total snow storage)
  snowpack$snow_total <- snowpack$snow_component + snowpack$liquid_component
  snowpack$snow_diff <- c(NA, -diff(snowpack$snow_total))
  snowpack$snow_diff[snowpack$snow_diff < 0] <- 0 
  
  #compute 7 day snowfall 
  snowpack$snow_sum_7days <- roll_sum(snowpack$snow_diff, n = 7, align = "right", fill = NA)
  
  #build precipitation time series and compute 7 day sum
  Rain_df <- data.frame(date = input_data$date, prec = input_data$prec)
  Rain_df <- slice(Rain_df, (warmup+1):nrow(Rain_df))
  Rain_df$rain_sum_7days <- roll_sum(Rain_df$prec, n = 7, align = "right", fill = NA)
  
  #build soil moisture storage time series and remove warmup period
  if (length_storage_df == n_dates){
    Soil_moisture_df <- data.frame(date=input_data$date, soil_moisture = storage_alle[,3])
  } else {
    Soil_moisture_df <- data.frame(date=input_data$date, soil_moisture = storage_alle[2:length_storage_df,3])
  }
  Soil_moisture_df <- slice(Soil_moisture_df, (warmup + 1):nrow(Soil_moisture_df))
  
  #indicator s7 snowmelt sum during the 7 day window ending on the flood date
  S7 <- c()
  for (d in flood_data_df$flood_date) {
    val <- snowpack[snowpack$date == as.Date(d), 6]
    S7 <- c(S7, val)
  }
  
  #indicator p7 precipitation sum during the 7 day window ending on the flood date
  P7 <- c()
  
  for (d in flood_data_df$flood_date) {
    val <- Rain_df[Rain_df$date == as.Date(d), 3]
    P7 <- c(P7, val)
  }
  
  #total water input indicator
  P_total <- S7 + P7
  
  #indicator pmax maximum daily precipitation within the 7 day window
  Pmax <- c()
  for (d in flood_data_df$flood_date) {
    window_start <- as.Date(d) - 6
    window_end <- as.Date(d)
    val <- max(Rain_df$prec[Rain_df$date >= window_start & Rain_df$date <= window_end], na.rm = TRUE)
    Pmax <- c(Pmax, val)
  }
  
  #indicator pmax frac fraction of the largest daily precipitation relative to total 7 day precipitation
  P_max_frac <- Pmax/P7
  
  #indicator sat start relative antecedent soil moisture at t0 (sm/fc)
  SM <- c()
  for (d in flood_data_df$t0) {
    val <- Soil_moisture_df[Soil_moisture_df$date == as.Date(d), 2]
    if (length(val)== 0){
      SM <- c (SM,NA)
    } else {
      SM <- c(SM, val)
    }
  }
  sat_start <- SM / as.numeric(param["FC"])   
  
  #compute thresholds used by the decision tree
  valid_snow <- snowpack$snow_sum_7days[snowpack$snow_sum_7days > 1]
  S7_quant <- quantile(valid_snow, 0.9, na.rm = TRUE)
  Pday_quant <- quantile(Rain_df$prec, 0.9, na.rm = TRUE) 
  P7_quant <- quantile(Rain_df$rain_sum_7days, 0.9, na.rm = TRUE)
  P7_mean <- mean(Rain_df$rain_sum_7days, na.rm = TRUE)

  #final output decision table 
  decision_df <-  data.frame(
    flood_date = flood_data_df$flood_date, 
    flood_magn = flood_data_df$flood_magn, 
    sat_start = round(sat_start, 2),
    Pmax = round(Pmax,2),
    P7 = round(P7,2),
    Pmax_frac = round(P_max_frac,2),
    S7 = round(S7,2),
    Ptotal = round(P_total,2),
    Pday_quant = round(Pday_quant,2), 
    P7_mean = round(P7_mean,2),
    P7_quant = round(P7_quant,2),
    S7_quant = round(S7_quant,2)
  )
  
  return(decision_df)
}


