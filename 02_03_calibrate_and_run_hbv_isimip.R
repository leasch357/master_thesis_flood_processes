##02_03 Load ISIMIP forcing, calibrate and run HBV model; calculate PET 

#install and load packages 
library(dplyr)

#fuctions 
source("R_func/hbv_sim_SAFE_updated.R")
source("R_func/hbv_snow_objfun_SAFE_updated.R")
source("R_func/pet_priestley_taylor.R")
source("R_func/calibrate_hbv_isimip.R")

#path to isimip forcing 
ordner_isimip <- "Output/ISIMIP/ISIMIP_forcing_hist" 

#output folder
ordner_isimip_output <- "Output/ISIMIP/ISIMIP_cal/" 

EZG <- read.csv2("Data/Camels/camels-20250724T0232Z/camels_name.txt", header = TRUE, colClasses = c("character", "character", "character")) #Ids EZG[,1]

#path to the CAMELS input data 
ordner_camels <- "Data/Camels/camels-20250724T0232Z/basin_timeseries_v1p2_modelOutput_daymet/model_output_daymet/model_output/flow_timeseries/daymet/" #01 muss variiert werden 
#load gauge id from camels 
EZG <- read.csv2("Data/Camels/camels-20250724T0232Z/camels_name.txt", header = TRUE, colClasses = c("character", "character", "character")) #Ids EZG[,1]

out_list <- lapply(
  X   = 1:nrow(EZG),
  FUN = calibrate_hbv_isimip,
  EZG = EZG,
  ordner_camels = ordner_camels,
  ordner_isimip = ordner_isimip,
  ordner_isimip_output = ordner_isimip_output
)
