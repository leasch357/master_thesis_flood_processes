#06 ISIMIP data evaluation: check climate bias (T, P, PET) by comparing CAMELS and ISIMIP forcing data 

#load packages 
library(dplyr)
library(lubridate)

#functions
source("R_func/pet_priestley_taylor.R")

#load catchment ids
EZG <- read.csv2("Data/CAMELS/camels-20250724T0232Z/camels_name.txt", header = TRUE, colClasses = c("character", "character", "character"))

#path to isimip forcing csv
ordner_isimip <- "Output/ISIMIP/"

#path to camels data 
ordner_camels <- "Data/CAMELS/camels-20250724T0232Z/basin_timeseries_v1p2_modelOutput_daymet/model_output_daymet/model_output/flow_timeseries/daymet/"

#output dataframe
clim_dev_df <- data.frame()

for (i in 1:nrow(EZG)) {
  g_id <- EZG [i,2]
  id <- EZG [i,1]
  
  #load camels time series for this catchment  
  pfad_cam <- paste0(g_id,"/",id,"_05_model_output.txt")
  data_cam <- read.table(paste0(ordner_camels, pfad_cam), header = TRUE)
  
  #create date column and extract data range
  data_cam$YR <- as.numeric(data_cam$YR)
  data_cam$MNTH <- as.numeric(data_cam$MNTH)
  data_cam$date <- as.Date(paste(data_cam$YR, data_cam$MNTH, data_cam$DY, sep = "-"))
  dateRange <- range(data_cam$date)
  
  #create camels data frame
  Camels_data <- data_cam[, c("date","PRCP","PET","OBS_RUN","TAIR")]
  names(Camels_data) <- c("date", "prec", "ept", "flow", "temp")
  
  #load historical isimip forcing data for this catchment
  Isimip_data_hist_all <- read.csv(paste0(ordner_isimip, "ISIMIP_forcing_hist/Forcing_df_", id, ".csv"))
  Isimip_data_hist_all$date <- as.Date(Isimip_data_hist_all$date)
  
  #adjust isimip to available camels time frame 
  Isimip_data_hist <- subset(Isimip_data_hist_all,
                             date >= dateRange[1] & date <= dateRange[2])
  
  #load parameters and select paramters with best run
  param_filepath <- paste0(ordner_isimip, "/ISIMIP_cal/all_runs_", id, "_ISIMIP.csv")
  param_isimip_all_runs <- read.csv(param_filepath)
  max_nse_index <- which.max(param_isimip_all_runs$NSE)
  param <- setNames(as.numeric(param_isimip_all_runs[max_nse_index, 1:13]), c("Ts","CFMAX","CFR","CWH","BETA","LP","FC","PERC","K0","K1","K2","UZL","MAXBAS"))
  alpha <- as.numeric(param_isimip_all_runs[max_nse_index, 14])
  
  #calculate PET 
  PET <- pet_priestley_taylor(forcing_df = Isimip_data_hist, alpha = alpha)
  
  Isimip_data <- data.frame(
    date = Isimip_data_hist$date,
    prec = Isimip_data_hist$prec,
    ept  = PET$PET,
    temp = Isimip_data_hist$temp
  )
  
  #temperature
  #mean temperature difference CAMELS and ISIMIP
  dT <- Isimip_data$temp - Camels_data$temp
  mean_dT <- mean(dT, na.rm = TRUE)
  mae_dT  <- mean(abs(dT), na.rm = TRUE)
  
  #mean monthly temperature difference 
  temp <- data.frame(
    date = Camels_data$date,
    month = month(Camels_data$date),
    temp_cam = Camels_data$temp,
    temp_isi = Isimip_data$temp
  )
  tmp_monthly <- temp %>%
    group_by(month) %>%
    summarise(
      id = id,
      temp_cam_m = mean(temp_cam, na.rm = TRUE),
      temp_isi_m = mean(temp_isi, na.rm = TRUE),
      temp_bias  = temp_isi_m - temp_cam_m,
      .groups = "drop"
    )
  
  mean_monthly_temp <- mean(tmp_monthly$temp_bias, na.rm = TRUE)
  mae_monthly_temp  <- mean(abs(tmp_monthly$temp_bias), na.rm = TRUE)
  
  #precipitation
  #mean precipitation difference
  dP <- Isimip_data$prec - Camels_data$prec
  mean_dP_mm <- mean(dP, na.rm = TRUE)
  mae_dP_mm  <- mean(abs(dP), na.rm = TRUE)
  
  #cumulative precipitation bias 
  sum_P_cam <- sum(Camels_data$prec, na.rm = TRUE)
  sum_P_isi <- sum(Isimip_data$prec, na.rm = TRUE)
  Pbias_pct <- 100 * (sum_P_isi - sum_P_cam) / sum_P_cam
  
  #mean annual cumulative precipitation bias 
  prec <- data.frame(
    year = year(Camels_data$date),
    P_cam = Camels_data$prec,
    P_isi = Isimip_data$prec
  )
    
  annual_bias <- prec %>%
    group_by(year) %>%
    summarise(
      sum_cam = sum(P_cam, na.rm=TRUE),
      sum_isi = sum(P_isi, na.rm=TRUE),
      bias_pct = 100 * (sum_isi - sum_cam) / sum_cam,
      .groups="drop"
    )
  
  mean_annual_bias_pct <- mean(annual_bias$bias_pct, na.rm = TRUE)
  mae_Pbias_annual_pct  <- mean(abs(annual_bias$bias_pct), na.rm = TRUE)
  
  #PET 
  #mean difference 
  dPET <- Isimip_data$ept - Camels_data$ept
  
  mean_dPET_mm <- mean(dPET, na.rm = TRUE)
  mae_dPET_mm  <- mean(abs(dPET), na.rm = TRUE)
  
  #mean monthly PET difference 
  pet <- data.frame(
    date = Camels_data$date,
    month = month(Camels_data$date),
    pet_cam = Camels_data$ept,
    pet_isi = Isimip_data$ept
  )
  
  pet_monthly <- pet %>%
    group_by(month) %>%
    summarise(
      id = id,
      pet_cam_m = mean(pet_cam, na.rm = TRUE),
      pet_isi_m = mean(pet_isi, na.rm = TRUE),
      pet_bias  = pet_isi_m - pet_cam_m,
      .groups = "drop"
    )
  
  mean_monthly_pet <- mean(pet_monthly$pet_bias, na.rm = TRUE)
  mae_monthly_pet  <- mean(abs(pet_monthly$pet_bias), na.rm = TRUE)
  
  #store in dataframe 
  clim_dev_df <- rbind(
    clim_dev_df,
    data.frame(
      id = id,
      mean_dT = mean_dT,
      mae_dT  = mae_dT,
      mean_monthly_temp = mean_monthly_temp, 
      mae_monthly_temp = mae_monthly_temp,
      
      mean_dP = mean_dP_mm,
      mae_dP  = mae_dP_mm,
      sum_P_cam = sum_P_cam,
      sum_P_isi = sum_P_isi,
      Pbias_cum_pct    = Pbias_pct,
      mean_annual_bias_pct = mean_annual_bias_pct,
      mae_Pbias_annual_pct = mae_Pbias_annual_pct,
      
      mean_dPET = mean_dPET_mm,
      mae_dPET  = mae_dPET_mm,
      mean_monthly_pet = mean_monthly_pet, 
      mae_monthly_pet = mae_monthly_pet
    )
  )
  
}

write.csv(clim_dev_df, "Output/ISIMIP/ISIMIP_output/clim_dev_df.csv")

