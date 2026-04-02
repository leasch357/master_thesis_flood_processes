#05 ISIMIP data evaluation: check mean occurence date Amax floods and magnitudes using circular statistics

library(dplyr)
library(lubridate)

source("R_func/hbv_sim_SAFE_updated.R")
source("R_func/hbv_snow_objfun_SAFE_updated.R")
source("R_func/pet_priestley_taylor.R")
source("R_Func/amax.R")
source("R_func/seasonality_functions.R")

#load catchment ids
EZG <- read.csv2("Data/CAMELS/camels-20250724T0232Z/camels_name.txt", header = TRUE, colClasses = c("character", "character", "character"))

#path to isimip forcing csv
ordner_isimip <- "Output/ISIMIP/"

#path to camels data 
ordner_camels <- "Data/CAMELS/camels-20250724T0232Z/basin_timeseries_v1p2_modelOutput_daymet/model_output_daymet/model_output/flow_timeseries/daymet/"

#define output dataframes 
Amax_df <- data.frame()
cir_stat_sum <- data.frame()
mag_sum <- data.frame()
var_80_00_hydro <- data.frame()

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
  param_filepath <- paste0(ordner_isimip, "ISIMIP_cal/all_runs_", id, "_ISIMIP.csv")
  param_isimip_all_runs <- read.csv(param_filepath)
  max_nse_index <- which.max(param_isimip_all_runs$NSE)
  param <- setNames(as.numeric(param_isimip_all_runs[max_nse_index, 1:13]), c("Ts","CFMAX","CFR","CWH","BETA","LP","FC","PERC","K0","K1","K2","UZL","MAXBAS"))
  alpha <- as.numeric(param_isimip_all_runs[max_nse_index, 14])
  
  #define warmup, calculate PET for warump period and built warmup df 
  warmup <- round(nrow(Camels_data)*0.1,0)
  start_main <- min(Camels_data$date)
  start_warm <- start_main - warmup
  
  Isimip_hist_warm <- subset(Isimip_data_hist_all, date >= start_warm & date < start_main)
  
  PET_warm <- pet_priestley_taylor(forcing_df = Isimip_hist_warm, alpha = alpha)
  
  Modelldata_warm <- data.frame(
    date = Isimip_hist_warm$date,
    prec = Isimip_hist_warm$prec,
    ept  = PET_warm$PET,
    flow = Camels_data$flow[1],
    temp = Isimip_hist_warm$temp
  )
  
  #calculate PET for main period and built modeldata df
  PET <- pet_priestley_taylor(forcing_df = Isimip_data_hist, alpha = alpha)
  
  Modelldata <- data.frame(
    date = Isimip_data_hist$date,
    prec = Isimip_data_hist$prec,
    ept  = PET$PET,
    flow = Camels_data$flow,
    temp = Isimip_data_hist$temp
  )
  
  #add warumup to main period
  Modelldata_isi_hist_inkl_warm <- rbind(Modelldata_warm, Modelldata)
  
  #run model
  model_output_isi_hist <- hbv_snow_objfun(
    param  = param,
    dat    = Modelldata_isi_hist_inkl_warm,
    warmup = warmup,
    Case   = 1
  )
  
  #Seasonality
  #Amax obs kalenderjahr
  Amax_obs <- Camels_data %>%
    mutate(jahr= year(date)) %>%
    group_by(jahr) %>% 
    slice_max(order_by = flow, n = 1, with_ties = FALSE) %>% 
    ungroup()%>%
    transmute(flood_date = as.Date(date), flood_magn = as.numeric(flow))
  
  Amax_isimip <- amax(input_data=Modelldata_isi_hist_inkl_warm, model_output=model_output_isi_hist,warmup=warmup)
  
  #erstes Jahr 1980 ist nur Oktober - Dez: hier kein Hochwasserevent für kal. jahr
  if (year(Amax_obs$flood_date[1]) == "1980"){
    Amax_obs    <- as.data.frame(Amax_obs[2:nrow(Amax_obs),])
    Amax_isimip <- as.data.frame(Amax_isimip[2:nrow(Amax_isimip),])
  }
  
  #Spearman correlation amax events per month 
  f_obs_amax <- get_rel_months(Amax_obs, date_col ="flood_date")
  f_isi_amax <- get_rel_months(Amax_isimip, date_col ="flood_date")
  r_amax_obs_isi_sp  <- cor(f_obs_amax, f_isi_amax, method = "spearman")
  
  Amax_df_id <- data.frame(
      id  = id,
      r_amax_obs_isi_sp = r_amax_obs_isi_sp
    )

  Amax_df <- rbind(Amax_df, Amax_df_id)

  #Circular statistics
  cir_stat_obs <- circ_stats_doy(date_col = Amax_obs$flood_date)
  cir_stat_isi <- circ_stats_doy(date_col = Amax_isimip$flood_date)
  
  mean_diff_days <- circ_diff_days(d1= cir_stat_obs$mean_doy, d2 = cir_stat_isi$mean_doy)

  cir_stat_sum_id <- data.frame(
    id                = id,
    mean_day_obs = cir_stat_obs$mean_day,
    mean_day_isi = cir_stat_isi$mean_day,
    diff_days    = mean_diff_days,
    R_obs = cir_stat_obs$R,
    R_isi= cir_stat_isi$R
  )

  cir_stat_sum <- rbind(cir_stat_sum, cir_stat_sum_id)

  #Magnitude
  mag_sum_id <- data.frame(
    id              = id,
    mean_mag_obs    = mean(Amax_obs$flood_magn),
    mean_mag_isi    = mean(Amax_isimip$flood_magn)
  )

  mag_sum_id$diff_mag <- ((mag_sum_id$mean_mag_isi - mag_sum_id$mean_mag_obs)/ mag_sum_id$mean_mag_obs) * 100
  
  mag_sum <- rbind(mag_sum, mag_sum_id)
}

write.csv(Amax_df, "Output/ISIMIP/ISIMIP_output/Amax_df.csv")
write.csv(cir_stat_sum, "Output/ISIMIP/ISIMIP_output/cir_stat_sum.csv")
write.csv(mag_sum, "Output/ISIMIP/ISIMIP_output/mag_sum.csv")
