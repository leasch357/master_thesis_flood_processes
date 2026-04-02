#03_02_Flood classification camels with both soil moisture calculation methods 

#load packages
library(lubridate)
library(dplyr)
library(RcppRoll)
library(changepoint)

#source founctions 
source("R_func/hbv_sim_SAFE_updated.R")
source("R_func/hbv_snow_objfun_SAFE_updated.R")
source("R_func/calculate_changepoint.R")
source("R_func/event_classification.R")
source("R_func/event_classification_CP_method.R")
source("R_func/event_classification_decision_tree.R")

#load catchment ids
EZG <- read.csv2("Data/CAMELS/camels-20250724T0232Z/camels_name.txt", header = TRUE, colClasses = c("character", "character", "character")) #Ids EZG[,1]

#path to camels topografy data
path_topo <- "Data/CAMELS/camels-20250724T0232Z/camels_topo.txt"

#path to inputdata
ordner <- "Data/CAMELS/camels-20250724T0232Z/basin_timeseries_v1p2_modelOutput_daymet/model_output_daymet/model_output/flow_timeseries/daymet/" ##01 muss variiert werden 

#load best parameter sets from camels calibration
param_camels <- read.csv("Output/CAMELS/Camels_cal/best_params_all_EZG_CAMELS.csv")

#define process labels and output dataframes 
all_processes <- c("rainandsnow", "snowmelt", "soilsat", "rainfall", "longrainfall", "noclass")
class_results_M1 <- data.frame()
class_results_M2 <- data.frame()

##Classification with fixed soil moisture threshold
#loop over all catchments
for (i in 1:nrow(EZG)) {
  
  g_id <- EZG [i,2]
  id <- EZG [i,1]
  
  #load camels time series for this catchment
  pfad <- paste0(g_id,"/",id,"_05_model_output.txt")
  Data <- read.table(paste0(ordner, pfad), header = TRUE)
  
  #create date column
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
  
  #prepare model input (forcing plus observed flow)
  Modelldata <- Data[time, c("date","PRCP","PET","OBS_RUN","TAIR")]
  names(Modelldata) <- c("date", "prec", "ept", "flow", "temp")
  
  #define warumup and add to input
  warmup <- round(nrow(Modelldata)*0.1,0)
  Modelldata_inkl_warm <- rbind(Modelldata [1:warmup,], Modelldata)
  
  Case <- 1 #Case = 1: interflow is dominant / Case = 2: percolation is dominant
  
  #load best parameter for this catchment
  param_opt <- param_camels[i,3:15]
  param_opt <- setNames(as.numeric(param_opt),c("Ts","CFMAX","CFR","CWH","BETA","LP","FC","PERC","K0","K1","K2","UZL","MAXBAS"))
  
  #run model
  model_output <- hbv_snow_objfun(param = param_opt, dat = Modelldata_inkl_warm, warmup = warmup, Case = Case)

  #create decision dataframe for decision tree with decision variables
  decision_df <- event_classification(
    input_data = Modelldata_inkl_warm,
    model_output = model_output,
    warmup = warmup,
    param = param_opt
  )
  
  #assign flood generating processes
  classification_result <- event_classification_quantile(decision_df = decision_df)
  
  #count process frequencies
  counts_M1 <- table(classification_result$mech_out)
  
  #ensure and add all processes are present in table 
  for (proc in all_processes) {
    if (!(proc %in% names(counts_M1))) counts_M1[proc] <- 0
  }
  counts_M1 <- counts_M1[all_processes]  # feste Reihenfolge
  
  #convert to dataframe
  counts_M1_df <- as.data.frame(as.list(counts_M1))
  
  #compute percentage shares
  perc_M1 <- round(100 * counts_M1_df[1, ] / sum(counts_M1_df[1, ]), 2)
  names(perc_M1) <- paste0("perc_", names(counts_M1_df))
  
  #identify dominant process and share
  dom_proc_M1  <- names(counts_M1_df)[which.max(counts_M1_df[1, ])]
  dom_share_M1 <- max(perc_M1)
  
  #store result for this catchment
  res_M1 <- data.frame(
    id = id,
    total_events = sum(counts_M1_df[1, ]),
    counts_M1_df,
    perc_M1,
    dominant_process = dom_proc_M1,
    dominant_share   = dom_share_M1
  )
  class_results_M1 <- rbind(class_results_M1, res_M1)
}

#safe results
write.csv(class_results_M1, "Output/CAMELS/CAMELS_output/classification_results_camels_M1.csv", row.names = FALSE)

##classification with changepoint-based soil moisture threshold 
#import camels changepoints and compute median changepoint for NA values
CP_df <- read.csv("Output/CAMELS/CAMELS_output/CP_CAMELS_full.csv")
CP_med <- median(CP_df$CP, na.rm = TRUE)
CP_df$id_char <- sprintf("%08d", CP_df$id)

#loop over all catchments
for (i in 1:nrow(EZG)) {
  
  g_id <- EZG [i,2]
  id <- EZG [i,1]
  
  #load camels time series for this catchment
  pfad <- paste0(g_id,"/",id,"_05_model_output.txt")
  Data <- read.table(paste0(ordner, pfad), header = TRUE)
  
  #create date column
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
  
  #prepare model input (forcing plus observed flow)
  Modelldata <- Data[time, c("date","PRCP","PET","OBS_RUN","TAIR")]
  names(Modelldata) <- c("date", "prec", "ept", "flow", "temp")
  
  #define warumup and add to input
  warmup <- round(nrow(Modelldata)*0.1,0)
  Modelldata_inkl_warm <- rbind(Modelldata [1:warmup,], Modelldata)
  
  Case <- 1
  
  #load best parameter for this catchment
  param_opt <- param_camels[i,3:15]
  param_opt <- setNames(as.numeric(param_opt),c("Ts","CFMAX","CFR","CWH","BETA","LP","FC","PERC","K0","K1","K2","UZL","MAXBAS"))
  
  #run model
  model_output <- hbv_snow_objfun(param = param_opt, dat = Modelldata_inkl_warm, warmup = warmup, Case = Case)
  
  #define CP for catchment
  CP_i <- CP_df$CP[CP_df$id_char == id]
  
  #run CP function for runoff peaks 
  threshold_SM <- calculate_changepoint(topo_path = path_topo, input_data = Modelldata_inkl_warm, model_output = model_output, warmup = warmup, id = id, show_plot = FALSE)
  
  #create decision dataframe for decision tree with decision variables
  decision_df_M2 <- event_classification_CP_method(input_data = Modelldata_inkl_warm, model_output = model_output, warmup = warmup, runoff_peaks = threshold_SM$runoff_peaks)
  
  #assign flood generating processes; if CP was NA the median is applied
  if (!is.na(CP_i)) {
    classification_result_M2 <- event_classification_quantile(
      decision_df = decision_df_M2, sat_thresh = CP_i)
    
  } else {
    classification_result_M2 <- event_classification_quantile(
      decision_df = decision_df_M2, sat_thresh = CP_med)
  }
  
  #count process frequencies
  counts_M2 <- table(classification_result_M2$mech_out)

  #ensure that all processes are present in table 
  for (proc in all_processes) {
    if (!(proc %in% names(counts_M2))) {
      counts_M2[proc] <- 0
    }
  }

  #convert to dataframe
  counts_M2_df <- as.data.frame(as.list(counts_M2))

  #compute percentage shares
  perc_M2 <- round(100 * counts_M2_df[1, ] / sum(counts_M2_df[1, ]), 2)
  names(perc_M2) <- paste0("perc_", names(counts_M2_df))
  
  #identify dominant process and share
  dom_proc_M2 <- names(counts_M2_df)[which.max(counts_M2_df[1, ])]
  dom_share_M2 <- max(perc_M2)

  #store result for this catchment
  res_M2 <- data.frame(id = id,
                       total_events = sum(counts_M2_df[1, ]),
                       counts_M2_df,
                       perc_M2,
                       dominant_process = dom_proc_M2,
                       dominant_share = dom_share_M2)
  class_results_M2 <- rbind(class_results_M2, res_M2)
}

#safe results
write.csv(class_results_M2, "Output/CAMELS/CAMELS_output/classification_results_camels_M2_test.csv", row.names = FALSE)






