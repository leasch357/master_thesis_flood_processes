calibrate_hbv_isimip <- function(i, EZG, ordner_camels, ordner_isimip, ordner_isimip_output){
  #calibrate HBV with ISIMIP input
  #Input:
  #i: row index of ezg
  #ezg: data frame with catchment id and gauge id
  #ordner_camels: path to camels basin time series files (used for observed runoff and date range)
  #ordner_isimip: path to isimip forcing csv files 
  #ordner_isimip_output: direction of output folder for calibration results
  #Output: 
  #list with best parameter row (max nse) and path to csv with all sampled runs

  g_id <- EZG [i,2]
  id <- EZG [i,1]
  
  #load camels data 
  pfad_cam <- paste0(g_id,"/",id,"_05_model_output.txt")
  data_cam <- read.table(paste0(ordner_camels, pfad_cam), header = TRUE)
  
  #define input data 
  data_cam$YR <- as.numeric(data_cam$YR)
  data_cam$MNTH <- as.numeric(data_cam$MNTH)
  data_cam$date <- as.Date(paste(data_cam$YR, data_cam$MNTH, data_cam$DY, sep = "-"))
  dateRange <- range(data_cam$date)
  Camels_data <- data_cam[, c("date","PRCP","PET","OBS_RUN","TAIR")]
  #rename column names 
  names(Camels_data) <- c("date", "prec", "ept", "flow", "temp")

  #load isimip forcing data  
  forcing_Isimip_data <- read.csv(paste0(ordner_isimip, "/Forcing_df_", id, ".csv"), stringsAsFactors = FALSE, check.names = FALSE)
  
  forcing_Isimip_data$date <- as.Date(forcing_Isimip_data$date)
  
  #clip isimip date to the same timeframe as camels
  Isimip_data <- subset(forcing_Isimip_data,
                       date >= dateRange[1] & date <= dateRange[2])
  
  ##define parametererange
  xmin <- c(-2, 0, 0, 0, 0.5, 0.3, 50, 0, 0.05, 0.01, 0.001, 0, 1, 1)
  xmax <- c(2, 10, 1, 0.2, 6, 1, 500, 5, 0.7, 0.5, 0.15, 100, 6, 1.6)
  
  param_names <- c("Ts","CFMAX", "CFR", "CWH", "BETA", "LP", "FC", "PERC", "K0", "K1", "K2","UZL", "MAXBAS", "alpha")
  
  #calibration runs
  n_sim <- 10000
  set.seed(1+i)  

  #dataframe with n_sim numbers of parameter combinations   
  param_samples <- data.frame(
    Ts     = runif(n_sim, min = xmin[1], max = xmax[1]),
    CFMAX  = runif(n_sim, min = xmin[2], max = xmax[2]),
    CFR    = runif(n_sim, min = xmin[3], max = xmax[3]),
    CWH    = runif(n_sim, min = xmin[4], max = xmax[4]),
    BETA   = runif(n_sim, min = xmin[5], max = xmax[5]),
    LP     = runif(n_sim, min = xmin[6], max = xmax[6]),
    FC     = runif(n_sim, min = xmin[7], max = xmax[7]),
    PERC   = runif(n_sim, min = xmin[8], max = xmax[8]),
    K0     = runif(n_sim, min = xmin[9], max = xmax[9]),
    K1     = runif(n_sim, min = xmin[10], max = xmax[10]),
    K2     = runif(n_sim, min = xmin[11], max = xmax[11]),
    UZL    = runif(n_sim, min = xmin[12], max = xmax[12]),
    MAXBAS = runif(n_sim, min = xmin[13], max = xmax[13]),
    alpha = runif(n_sim, min = xmin[14], max = xmax[14])
  )
  
  colnames(param_samples) <- param_names
  param_samples$NSE <- NA_real_
  
  #calibration of HBV model 
  
    for (j in 1:n_sim) {
      param <- as.numeric(param_samples[j,1:13])
      alpha <- as.numeric(param_samples[j,14])
      
      #warmup
      warmup_years <- 3
      start_warm <- dateRange[1] - 365*warmup_years
      end_warm <- dateRange[1]-1
      isimip_warm <- subset(forcing_Isimip_data, date >= start_warm & date <= end_warm)
      
      PET_warm <- pet_priestley_taylor(forcing_df = isimip_warm, alpha = alpha)
      
      Modelldata_isimip_warm <- data.frame(
        date = isimip_warm$date,
        prec = isimip_warm$prec,
        ept  = PET_warm$PET,
        flow = Camels_data$flow[1],
        temp = isimip_warm$temp
      )
      
      PET <- pet_priestley_taylor(forcing_df = Isimip_data, alpha = alpha)
      
      Modelldata_isimip <- data.frame(
        date = Isimip_data$date,
        prec = Isimip_data$prec,
        ept  = PET$PET,
        flow = Camels_data$flow,     
        temp = Isimip_data$temp
      )
      
      #warumup period is placed at beginning of model data 
      Modelldata_hbv_isimip <- rbind(Modelldata_isimip_warm, Modelldata_isimip)
      warmup_days <- nrow(Modelldata_isimip_warm)
      
      model_output <- hbv_snow_objfun(param = param, dat = Modelldata_hbv_isimip, warmup = warmup_days, Case = 1)
      
      nse <- model_output[2]
      
      param_samples$NSE[j] <- nse
    }
  
  out_file <- paste0(ordner_isimip_output, "all_runs_", id, "_ISIMIP.csv")
  write.csv(param_samples, out_file, row.names = FALSE)
  
  res <- list(
    out_file = out_file
  )
  
  return(res)
}
  
   