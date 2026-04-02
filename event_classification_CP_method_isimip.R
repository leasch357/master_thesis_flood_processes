event_classification_CP_method_isimip <- function(input_data, model_output,warmup, runoff_peaks, start_date = NULL, end_date = NULL) {
  #prepare decision variables for flood process classification
  #method 2 uses normalized antecedent soil moisture based on runoff peaks (changepoint method)
  #input:
  #input_data model forcing and observed flow
  #model_output hbv output containing qsim and states
  #warmup number of warmup days to remove
  #runoff_peaks include sm values used for min max normalization
  #start_date and end_date is used for classification of specific periods e.g. future periods 
  #output:
  #decision df with flood variables (indicators, thresholds)
  
  #extract simulated discharge and remove warmup period
  Q_sim <- attr(model_output, "Q_sim")
  output_data <- data.frame(date = input_data$date, Q_sim = Q_sim)
  output_data <- slice(output_data, (warmup + 1):nrow(output_data))
  
  #ensure series starts on jan 1 by dropping the first incomplete year
  first_year <- year(output_data$date[1])
  if (yday(output_data$date[1]) != 1) {
    output_data <- output_data %>% filter(year(date) > first_year)}
  
  #extract annual maxima (amax) of simulated discharge
  flood_data_df <- output_data %>%
    mutate(jahr= year(date)) %>%
    group_by(jahr) %>% 
    slice_max(order_by = Q_sim, n = 1) %>% 
    ungroup()
  
  #keep date and magnitude and define t0 (start of 7 day window)
  flood_data_df <- flood_data_df[, c("date", "Q_sim")]
  colnames(flood_data_df) = c("flood_date", "flood_magn")
  flood_data_df$t0 <- flood_data_df$flood_date - 6 
  
  #adjust to selected period
  if (!is.null(start_date)) {
    flood_data_df <- flood_data_df %>% filter(flood_date >= start_date)
  }
  if (!is.null(end_date)) {
    flood_data_df <- flood_data_df %>% filter(flood_date <= end_date)
  }
  
  #extract hbv state variables (snow and soil moisture) 
  storage_alle <- attr(model_output, "STATES")
  length_storage_df <- length(storage_alle[,1])
  n_dates  <- nrow(input_data)
  
  #build snow storages time series and remove warmup period
  if (length_storage_df == n_dates){
    snowpack <- data.frame(date=input_data$date, snow_component = storage_alle[,1], liquid_component = storage_alle[,2])
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
  
  #adjust to selected period
  if (!is.null(start_date)) {
    snowpack <- snowpack %>% dplyr::filter(date >= start_date)
  }
  if (!is.null(end_date)) {
    snowpack <- snowpack %>% dplyr::filter(date <= end_date)
  }
  
  #build precipitation time series and compute 7 day sum
  Rain_df <- data.frame(date = input_data$date, prec = input_data$prec)
  Rain_df <- slice(Rain_df, (warmup+1):nrow(Rain_df))
  Rain_df$rain_sum_7days <- roll_sum(Rain_df$prec, n = 7, align = "right", fill = NA)
  
  #adjust to selected period
  if (!is.null(start_date)) {
    Rain_df <- Rain_df %>% dplyr::filter(date >= start_date-7)
  }
  if (!is.null(end_date)) {
    Rain_df <- Rain_df %>% dplyr::filter(date <= end_date)
  }
  
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
  
  #indicator pmax maximum daily precipitation within the 7 day window
  P7 <- c()
  for (d in flood_data_df$flood_date) {
    val <- Rain_df[Rain_df$date == as.Date(d), 3]
    P7 <- c(P7, val)
  }
  
  #total water input indicator
  P_total <- S7 + P7
  
  #indicator pmax frac fraction of the largest daily precipitation relative to total 7 day precipitation
  Pmax <- c()
  for (d in flood_data_df$flood_date) {
    window_start <- as.Date(d) - 6
    window_end <- as.Date(d)
    val <- max(Rain_df$prec[Rain_df$date >= window_start & Rain_df$date <= window_end], na.rm = TRUE)
    Pmax <- c(Pmax, val)
  }
  
  #indicator pmax frac fraction of the largest daily precipitation relative to total 7 day precipitation
  P_max_frac <- Pmax/P7
  
  #indicator sat start relative antecedent soil moisture at t0
  SM <- c()
  for (d in flood_data_df$t0) {
    val <- Soil_moisture_df[Soil_moisture_df$date == as.Date(d), 2]
    if (length(val)== 0){
      SM <- c (SM,NA)
    } else {
      SM <- c(SM, val)
    }
  }
  
  #normalize soil moisture values
  min_SM <- min(runoff_peaks$SM, na.rm = TRUE)
  max_SM <- max(runoff_peaks$SM, na.rm = TRUE)
  sat_start_norm <- (SM - min_SM) / (max_SM - min_SM)
  
  #compute thresholds used by the decision tree
  valid_snow <- snowpack$snow_sum_7days[snowpack$snow_sum_7days > 1]
  S7_quant <- as.numeric(quantile(valid_snow, 0.9, na.rm = TRUE))
  Pday_quant <- as.numeric(quantile(Rain_df$prec, 0.9, na.rm = TRUE))
  P7_quant <- as.numeric(quantile(Rain_df$rain_sum_7days, 0.9, na.rm = TRUE))
  P7_mean <- mean(Rain_df$rain_sum_7days, na.rm = TRUE) 
  
  #final output decision table  
  decision_df <-  data.frame(
    flood_date = flood_data_df$flood_date, 
    flood_magn = flood_data_df$flood_magn, 
    sat_start = round(sat_start_norm, 2),
    Pmax = round(Pmax,2),
    P7 = round(P7,2),
    Pmax_frac = round(P_max_frac,2),
    S7 = round(S7,2),
    Ptotal = round(P_total,2),
    Pday_quant = round(Pday_quant,2), #ab hier Schwellenwerte
    P7_mean = round(P7_mean,2),
    P7_quant = round(P7_quant,2),
    S7_quant = round(S7_quant,2)
  )
  
  return(decision_df)
}

