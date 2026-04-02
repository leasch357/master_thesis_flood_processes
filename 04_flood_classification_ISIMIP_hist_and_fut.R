#04_Flood classification isimip for historical and future period

#load packages 
library(dplyr)
library(RcppRoll)
library(lubridate)
library(stats)
library(changepoint)

#source founctions 
source("R_func/hbv_sim_SAFE_updated.R")
source("R_func/hbv_snow_objfun_SAFE_updated.R")
source("R_func/pet_priestley_taylor.R")
source("R_func/calculate_changepoint.R")
source("R_func/event_classification_CP_method_isimip.R")
source("R_func/event_classification_decision_tree.R")

#load catchment ids
EZG <- read.csv2("Data/CAMELS/camels-20250724T0232Z/camels_name.txt", header = TRUE, colClasses = c("character", "character", "character"))

#path to camels topografy data
path_topo <- "Data/CAMELS/camels-20250724T0232Z/camels_topo.txt"

#path to isimip forcing csv
ordner_isimip <- "Output/ISIMIP/"

#path to camels data 
ordner_camels <- "Data/CAMELS/camels-20250724T0232Z/basin_timeseries_v1p2_modelOutput_daymet/model_output_daymet/model_output/flow_timeseries/daymet/"

#load changepoints and compute median changepoint
CP_cam_hist <- read.csv("Output/CAMELS/CAMELS_output/CP_CAMELS_full.csv")
sm_med <- median(CP_cam_hist$CP, na.rm = TRUE)
CP_cam_hist$id_char <- sprintf("%08d", CP_cam_hist$id)

#define process labels and output dataframes 
all_processes <- c("rainandsnow", "snowmelt", "soilsat", "rainfall", "longrainfall", "noclass")
class_results_hist  <- data.frame()
class_results_future_41_70  <- data.frame()
class_results_future_71_00  <- data.frame() 
flood_proc_hist  <- data.frame()
flood_proc_future_41_70  <- data.frame()
flood_proc_future_71_00  <- data.frame()

#loop over all catchments
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
  
  #historical period
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
  
  #processes historical
  #initialize placeholders for dominant process summaries
  Dom_flood_proc_cam <- NA
  Prozent_cam        <- NA
  Dom_flood_proc_isi <- NA
  Prozent_isi        <- NA
  
  #select catchment specific soil moisture threshold from camels changepoints
  CP <- CP_cam_hist$CP[CP_cam_hist$id_char == id]

  #derive runoff peaks and soil moisture threshold for isimip historical simulation
  sm_thresh_isi <- calculate_changepoint(topo_path = path_topo,input_data=Modelldata_isi_hist_inkl_warm, model_output=model_output_isi_hist,warmup=warmup, id=id, show_plot = FALSE)
  
  #build decision data frame for decision tree using runoff peaks
  decision_df_isimip_method2 <- event_classification_CP_method_isimip(input_data=Modelldata_isi_hist_inkl_warm, model_output=model_output_isi_hist,warmup=warmup, runoff_peaks = sm_thresh_isi$runoff_peaks)
  
  #classify flood generating processes using catchment specific cp if available, otherwise median threshold
  if (!is.na(CP)) {
    
    classification_result_isimip_method2 <- event_classification_quantile(decision_df=decision_df_isimip_method2, parts_rainsnow = 1/3, parts_fracextreme = 2/3, sat_thresh = CP)
    
  } else {
    classification_result_isimip_method2 <- event_classification_quantile(decision_df=decision_df_isimip_method2, parts_rainsnow = 1/3, parts_fracextreme = 2/3, sat_thresh = sm_med)
  }
  
  #count process frequencies and ensure all process labels are present
  counts_M2 <- table(classification_result_isimip_method2$mech_out)
  for (proc in all_processes) {
    if (!(proc %in% names(counts_M2))) {
      counts_M2[proc] <- 0
    }
  }
  
  #convert counts to data frame and compute percentage shares
  counts_M2_df <- as.data.frame(as.list(counts_M2))
  perc_M2 <- round(100 * counts_M2_df[1, ] / sum(counts_M2_df[1, ]), 2)
  names(perc_M2) <- paste0("perc_", names(counts_M2_df))
  
  #identify dominant process and its share
  dom_proc_M2 <- names(counts_M2_df)[which.max(counts_M2_df[1, ])]
  dom_share_M2 <- max(perc_M2)
  
  #store full classification summary for this catchment (historical)
  res_M2 <- data.frame(id = id,
                       total_events = sum(counts_M2_df[1, ]),
                       counts_M2_df,
                       perc_M2,
                       dominant_process = dom_proc_M2,
                       dominant_share = dom_share_M2)
  class_results_hist <- rbind(class_results_hist, res_M2)
  
  #store only dominant process and share for compact output table
  dom_flood_class_isi <- classification_result_isimip_method2 %>%
    count(mech_out) %>%
    mutate(Prozent = round(100 * n / sum(n), 1)) %>%
    arrange(desc(n))%>%
    slice(1)
  
  Dom_flood_proc_isi <- dom_flood_class_isi[1,1]
  Prozent_isi        <- dom_flood_class_isi[1,3]
  
  flood_proc_future_id <- data.frame(
    id                 = id,
    Dom_flood_proc_isi = Dom_flood_proc_isi,
    Prozent_isi = Prozent_isi
  )
  
  flood_proc_hist <- rbind (flood_proc_hist,flood_proc_future_id)
  
  #future period
  #load isimip future forcing ssp370/585
  Isimip_data_fut <- read.csv(paste0(ordner_isimip, "/ISIMIP_forcing_future/ssp370/Forcing_df_", id, "_ssp370.csv"))
  #Isimip_data <- read.csv(paste0(ordner_isimip, "/ISIMIP_forcing_future/Forcing_df_", id, ".csv")) #ssp558
  Isimip_data_fut$date <- as.Date(Isimip_data_fut$date)
  
  #define warmup period for future simulation (here 3 years)
  warmup_years <- 3
  dateRange <- range(Isimip_data_fut$date)
  start_warm <- dateRange[1] - 365*warmup_years
  end_warm <- dateRange[1]-1
  isimip_warm <- subset(Isimip_data_hist_all, date >= start_warm & date <= end_warm)
  
  #calculate pet for warmup and build warmup model input
  PET_warm <- pet_priestley_taylor(forcing_df = isimip_warm, alpha = alpha)
  Modelldata_isimip_warm <- data.frame(
    date = isimip_warm$date,
    prec = isimip_warm$prec,
    ept  = PET_warm$PET,
    flow = rep(NA_real_, nrow(isimip_warm)),
    temp = isimip_warm$temp
  )
  
  #calculate pet for main future period and build model input
  PET <- pet_priestley_taylor(forcing_df = Isimip_data_fut, alpha = alpha)
  Modelldata_isimip <- data.frame(
    date = Isimip_data_fut$date,
    prec = Isimip_data_fut$prec,
    ept  = PET$PET,
    flow = rep(NA_real_, nrow(Isimip_data_fut)),
    temp = Isimip_data_fut$temp
  )
  
  #prepend warmup to future period
  Modelldata_isi_fut_inkl_warm <- rbind(Modelldata_isimip_warm, Modelldata_isimip)
  warmup_isi <- nrow(Modelldata_isimip_warm)
  
  #run hbv
  model_output_isimip_fut <- hbv_snow_objfun(param = param, dat = Modelldata_isi_fut_inkl_warm, warmup = warmup_isi, Case = 1)
  
  #initialize placeholders for dominant process summaries (future)
  Dom_flood_proc_cam <- NA
  Prozent_cam        <- NA
  Dom_flood_proc_isi <- NA
  Prozent_isi        <- NA
  
  #future period from 2041-2070
  #derive runoff peaks for future 
  sm_thresh_isi <- calculate_changepoint(topo_path = path_topo,input_data=Modelldata_isi_fut_inkl_warm, model_output=model_output_isimip_fut,warmup=warmup_isi, id=id, show_plot = FALSE)
  
  #build decision data frame for decision tree for selected future subperiod
  decision_df_isimip_method2 <- event_classification_CP_method_isimip(input_data=Modelldata_isi_fut_inkl_warm, model_output=model_output_isimip_fut,warmup=warmup_isi, runoff_peaks = sm_thresh_isi$runoff_peaks, start_date = as.Date("2041-01-01"), end_date  = as.Date("2070-12-31"))
  
  #classify flood generating processes using catchment specific cp if available, otherwise median threshold
  if (!is.na(CP)) {
    classification_result_isimip_method2 <- event_classification_quantile(decision_df=decision_df_isimip_method2, parts_rainsnow = 1/3, parts_fracextreme = 2/3, sat_thresh = CP)
    
  } else {
    classification_result_isimip_method2 <- event_classification_quantile(decision_df=decision_df_isimip_method2, parts_rainsnow = 1/3, parts_fracextreme = 2/3, sat_thresh = sm_med)
  }
  
  #summarize classification results for this catchment and subperiod
  counts_M2 <- table(classification_result_isimip_method2$mech_out)
  for (proc in all_processes) {
    if (!(proc %in% names(counts_M2))) {
      counts_M2[proc] <- 0
    }
  }
  counts_M2_df <- as.data.frame(as.list(counts_M2))
  perc_M2 <- round(100 * counts_M2_df[1, ] / sum(counts_M2_df[1, ]), 2)
  names(perc_M2) <- paste0("perc_", names(counts_M2_df))
  dom_proc_M2 <- names(counts_M2_df)[which.max(counts_M2_df[1, ])]
  dom_share_M2 <- max(perc_M2)
  
  res_M2 <- data.frame(id = id,
                       total_events = sum(counts_M2_df[1, ]),
                       counts_M2_df,
                       perc_M2,
                       dominant_process = dom_proc_M2,
                       dominant_share = dom_share_M2)
  
  class_results_future_41_70 <- rbind(class_results_future_41_70, res_M2)
  
  dom_flood_class_isi <- classification_result_isimip_method2 %>%
    count(mech_out) %>%
    mutate(Prozent = round(100 * n / sum(n), 1)) %>%
    arrange(desc(n))%>%
    slice(1)
  
  Dom_flood_proc_isi <- dom_flood_class_isi[1,1]
  Prozent_isi        <- dom_flood_class_isi[1,3]
  
  flood_proc_future_id <- data.frame(
    id                 = id,
    Dom_flood_proc_isi = Dom_flood_proc_isi,
    Prozent_isi = Prozent_isi
  )
  
  flood_proc_future_41_70 <- rbind (flood_proc_future_41_70,flood_proc_future_id)
  
  #future period from 2071-2100
  decision_df_isimip_method2 <- event_classification_CP_method_isimip(input_data=Modelldata_isi_fut_inkl_warm, model_output=model_output_isimip_fut,warmup=warmup_isi, runoff_peaks = sm_thresh_isi$runoff_peaks, start_date = as.Date("2071-01-01"), end_date  = as.Date("2100-12-31"))
  
  if (!is.na(CP)) {
    classification_result_isimip_method2 <- event_classification_quantile(decision_df=decision_df_isimip_method2, parts_rainsnow = 1/3, parts_fracextreme = 2/3, sat_thresh = CP)
    
  } else {
    classification_result_isimip_method2 <- event_classification_quantile(decision_df=decision_df_isimip_method2, parts_rainsnow = 1/3, parts_fracextreme = 2/3, sat_thresh = sm_med)
  }
  
  counts_M2 <- table(classification_result_isimip_method2$mech_out)
  for (proc in all_processes) {
    if (!(proc %in% names(counts_M2))) {
      counts_M2[proc] <- 0
    }
  }
  counts_M2_df <- as.data.frame(as.list(counts_M2))
  perc_M2 <- round(100 * counts_M2_df[1, ] / sum(counts_M2_df[1, ]), 2)
  names(perc_M2) <- paste0("perc_", names(counts_M2_df))
  dom_proc_M2 <- names(counts_M2_df)[which.max(counts_M2_df[1, ])]
  dom_share_M2 <- max(perc_M2)
  
  res_M2 <- data.frame(id = id,
                       total_events = sum(counts_M2_df[1, ]),
                       counts_M2_df,
                       perc_M2,
                       dominant_process = dom_proc_M2,
                       dominant_share = dom_share_M2)
  
  class_results_future_71_00 <- rbind(class_results_future_71_00, res_M2)
  
  dom_flood_class_isi <- classification_result_isimip_method2 %>%
    count(mech_out) %>%
    mutate(Prozent = round(100 * n / sum(n), 1)) %>%
    arrange(desc(n))%>%
    slice(1)
  
  Dom_flood_proc_isi <- dom_flood_class_isi[1,1]
  Prozent_isi        <- dom_flood_class_isi[1,3]
  
  flood_proc_future_id <- data.frame(
    id                 = id,
    Dom_flood_proc_isi = Dom_flood_proc_isi,
    Prozent_isi = Prozent_isi
  )
  
  flood_proc_future_71_00 <- rbind (flood_proc_future_71_00,flood_proc_future_id)
}

write.csv(flood_proc_hist, "Output/ISIMIP/ISIMIP_output/flood_proc_hist_CP_camels.csv")
write.csv(flood_proc_future_41_70, "Output/ISIMIP/ISIMIP_output/flood_proc_future_2041_2070_ssp370_CP_camels.csv")
write.csv(flood_proc_future_71_00, "Output/ISIMIP/ISIMIP_output/flood_proc_future_2071_2100_ssp370_CP_camels.csv")

write.csv(class_results_hist, "Output/ISIMIP/ISIMIP_output/class_results_hist_CP_camels.csv")
write.csv(class_results_future_41_70, "Output/ISIMIP/ISIMIP_output/class_results_future_2041_2070_ssp370_CP_camels.csv")
write.csv(class_results_future_71_00, "Output/ISIMIP/ISIMIP_output/class_results_future_2071_2100_ssp370_CP_camels.csv")

