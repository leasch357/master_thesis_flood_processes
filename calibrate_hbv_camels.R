calibrate_hbv_cam <- function(i, EZG, ordner_camels, ordner_camels_output){
  #calibrate HBV with CAMELS input
  #Input:
  #i: row index of ezg
  #ezg: data frame with catchment id and gauge id
  #ordner_camels: path to camels basin time series files 
  #ordner_camels_output: direction of output folder for calibration results
  #Output: 
  #list with best parameter row (max nse) and path to csv with all sampled runs

  g_id <- EZG [i,2]
  id <- EZG [i,1]
  
  #load CAMELS
  pfad_cam <- paste0(g_id,"/",id,"_05_model_output.txt")
  data_cam <- read.table(paste0(ordner_camels, pfad_cam), header = TRUE)
  
  #define input data 
  data_cam$YR <- as.numeric(data_cam$YR)
  data_cam$MNTH <- as.numeric(data_cam$MNTH)
  data_cam$date <- as.Date(paste(data_cam$YR, data_cam$MNTH, data_cam$DY, sep = "-"))
  dateRange <- range(data_cam$date)
  camels_data <- data_cam[, c("date","PRCP","PET","OBS_RUN","TAIR")]
  
  #rename column names 
  names(camels_data) <- c("date", "prec", "ept", "flow", "temp")
  
  #define warmup 
  warmup <- round(nrow(camels_data)*0.1,0)
  modelinput_camels <- rbind(camels_data [1:warmup,], camels_data)
  
  #define parametererange
  xmin <- c(-2, 0, 0, 0, 0.5, 0.3, 50, 0, 0.05, 0.01, 0.001, 0, 1)
  xmax <- c(2, 10, 1, 0.2, 6, 1, 500, 5, 0.7, 0.5, 0.15, 100, 6)
  
  param_names <- c("Ts","CFMAX", "CFR", "CWH", "BETA", "LP", "FC", "PERC", "K0", "K1", "K2","UZL", "MAXBAS")
  
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
    MAXBAS = runif(n_sim, min = xmin[13], max = xmax[13])
  )
  
  colnames(param_samples) <- param_names
  param_samples$NSE <- NA
  
  #calibration of HBV model 
  
  for (j in 1:n_sim) {
    param <- as.numeric(param_samples[j,1:13])

    model_output <- hbv_snow_objfun(param = param, dat = modelinput_camels, warmup = warmup, Case = 1) 
    
    nse <- model_output[2]
    
    param_samples$NSE[j] <- nse
  }
  
  out_file <- paste0(ordner_camels_output, "all_runs_", id, "_CAMELS.csv")
  write.csv(param_samples, out_file, row.names = FALSE)
  
  res <- list(
    out_file = out_file
  )
  
  return(res)
}

