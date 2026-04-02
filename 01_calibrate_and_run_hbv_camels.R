##01 Load CAMELS, calibrate and run HBV model

#set working directory
setwd("C:/Users/leasc/Desktop/Masterarbeit/R/R_Abgabe")

#install and load packages 
library(dplyr)

#functions
source("R_func/hbv_sim_SAFE_updated.R")
source("R_func/hbv_snow_objfun_SAFE_updated.R")
source("R_func/calibrate_hbv_camels.R")

#output folder
ordner_camels_output <- "Output/CAMELS/CAMELS_cal/" 

#path to the CAMELS input data 
ordner_camels <- "Data/CAMELS/camels-20250724T0232Z/basin_timeseries_v1p2_modelOutput_daymet/model_output_daymet/model_output/flow_timeseries/daymet/"

#load gauge id 
EZG <- read.csv2("Data/Camels/camels-20250724T0232Z/camels_name.txt", header = TRUE, colClasses = c("character", "character", "character")) #Ids EZG[,1]

out_list <- lapply(
  X   = 1:nrow(EZG),
  FUN = calibrate_hbv_cam,
  EZG = EZG,
  ordner_camels = ordner_camels,
  ordner_camels_output = ordner_camels_output
)

