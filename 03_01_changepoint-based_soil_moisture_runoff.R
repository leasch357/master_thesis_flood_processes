#03_01 Derive soil-moisture changepoints from CAMELS data 

#load required packages
library(dplyr)
library(changepoint)
library(ggplot2)

#source HBV functions and changepoint routine 
source("R_func/hbv_sim_SAFE_updated.R")
source("R_func/hbv_snow_objfun_SAFE_updated.R")
source("R_func/calculate_changepoint.R")

#path to Camels input data 
ordner_camels_input <- "Data/CAMELS/camels-20250724T0232Z/basin_timeseries_v1p2_modelOutput_daymet/model_output_daymet/model_output/flow_timeseries/daymet/"

#path to camels topografy data
path_topo <- "Data/CAMELS/camels-20250724T0232Z/camels_topo.txt"

#load catchment list with ids and gauge ids
EZG <- read.csv2("Data/CAMELS/camels-20250724T0232Z/camels_name.txt", header = TRUE, colClasses = c("character", "character", "character")) #Ids EZG[,1]

#load best parameter sets from camels calibration
param_camels <- read.csv("Output/CAMELS/Camels_cal/best_params_all_EZG_CAMELS.csv", header = TRUE)

#create empty output dataframe to store changepoint thresholds
CP_camels <- data.frame()

#directory to to save changpoint plots 
plot_dir <- "Output/CAMELS/Plots_changepoint"

#loop over all catchemnts
for (i in 1:nrow(EZG)) {
  
  g_id <- EZG [i,2]
  id <- EZG [i,1]
  
  #load CAMELS time series for this catchment
  pfad <- paste0(g_id,"/",id,"_05_model_output.txt")
  Data <- read.table(paste0(ordner_camels_input, pfad), header = TRUE)
  
  #create data column
  Data$YR <- as.numeric(Data$YR)
  Data$MNTH <- as.numeric(Data$MNTH)
  Data$HR <- as.numeric(Data$HR)
  Data$date <- as.Date(paste(Data$YR, Data$MNTH, Data$DY, sep = "-"))
  
  #define available time period
  min_date <- min(Data$date)
  max_date <- max(Data$date)
  dateRange <- c(min_date,max_date)
  
  t_start <- which(Data$date == dateRange[1])
  t_end <- which(Data$date == dateRange[2])
  time <- t_start:t_end
  
  #prepare model input data
  Modelldata <- Data[time, c("date","PRCP","PET","OBS_RUN","TAIR")]
  names(Modelldata) <- c("date", "prec", "ept", "flow", "temp")
  
  #define warmup-period (first 10% of record) and add it to the model input
  warmup <- round(nrow(Modelldata)*0.1,0)
  Modelldata_inkl_warm <- rbind(Modelldata [1:warmup,], Modelldata)
  
  Case <- 1 #Case = 1: interflow is dominant / Case = 2: percolation is dominant
  
  #load calibrated parameter set for this catchment
  param_opt <- param_camels[i,3:15]
  param_opt <- setNames(as.numeric(param_opt),c("Ts","CFMAX","CFR","CWH","BETA","LP","FC","PERC","K0","K1","K2","UZL","MAXBAS"))
  
  #run model
  model_output <- hbv_snow_objfun(param = param_opt, dat = Modelldata_inkl_warm, warmup = warmup, Case = Case)
  
  #changepoint-based soil moisture threshold 
  threshold_SM <- calculate_changepoint(topo_path = path_topo, input_data = Modelldata_inkl_warm, model_output = model_output, warmup = warmup, id = id, show_plot = TRUE)
  
  #skip catchment if no changepoint was detected
  if (isFALSE(threshold_SM$has_cp)) next
  
  #store changepoint threshold and model fit metric (nRMSE)
  CP <- data.frame(
    id = id, 
    CP = threshold_SM$threshold_SM,
    nRMSE = threshold_SM$nRMSE
  )
  
  CP_camels <- rbind(CP_camels,CP)
  
  #save changepoint plot
  if (!is.null(threshold_SM$plot)){
    ggsave(filename = paste0(plot_dir, "/", id, "_SM_threshold.png"),
           plot = threshold_SM$plot,
           width = 7, height = 5, dpi = 300,
           bg = "white")}
  
}

#save results
write.csv(CP_camels, "Output/CAMELS/CAMELS_output/CP_CAMELS_full.csv", row.names = FALSE)



