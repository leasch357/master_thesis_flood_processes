##02_01 Extraction of isimip raster data per catchment; creating of forcing data set for historical period 
#load packages
library(terra)

#functions
source("R_func/raster_extraction.R")

#output folder
ordner_isimip_output <- "Output/ISIMIP/ISIMIP_forcing_hist"
scenario <- "hist"

#path to historical ISIMIP raster data 
pfad <- "Data/ISIMIP_hist"

#load catchment shapes and catchment ids 
EZG_shape <- vect("Data/Camels/camels-20250724T0232Z/basin_set_full_res/HCDN_nhru_final_671.shp")
EZG <- read.csv2("Data/Camels/camels-20250724T0232Z/camels_name.txt", header = TRUE, colClasses = c("character", "character", "character"))
EZG_id <- EZG[,1]

#load raster data
#temperature
r_tas_71_80 <- rast(file.path(
  pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_tas_global_daily_1971_1980.nc"))
r_tas_81_90 <- rast(file.path(   
  pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_tas_global_daily_1981_1990.nc"))
r_tas_91_00 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_tas_global_daily_1991_2000.nc"))
r_tas_01_10 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_tas_global_daily_2001_2010.nc"))
r_tas_11_14 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_tas_global_daily_2011_2014.nc"))
r_tas <- c(r_tas_71_80, r_tas_81_90, r_tas_91_00, r_tas_01_10, r_tas_11_14)
time(r_tas) <- c(time(r_tas_71_80), time(r_tas_81_90), time(r_tas_91_00), time(r_tas_01_10), time(r_tas_11_14))

#precipitation
r_pr_71_80 <- rast(file.path(   
  pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_pr_global_daily_1971_1980.nc"))
r_pr_81_90 <- rast(file.path(   
  pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_pr_global_daily_1981_1990.nc"))
r_pr_91_00 <- rast(file.path(   
  pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_pr_global_daily_1991_2000.nc"))
r_pr_01_10 <- rast(file.path(   
  pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_pr_global_daily_2001_2010.nc"))
r_pr_11_14 <- rast(file.path(   
  pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_pr_global_daily_2011_2014.nc"))
r_pr <- c(r_pr_71_80, r_pr_81_90, r_pr_91_00, r_pr_01_10, r_pr_11_14)
time(r_pr) <- c(time(r_pr_71_80), time(r_pr_81_90), time(r_pr_91_00), time(r_pr_01_10), time(r_pr_11_14))

#air pressure
r_ps_71_80 <- rast(file.path(   
  pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_ps_global_daily_1971_1980.nc"))
r_ps_81_90 <- rast(file.path(   
  pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_ps_global_daily_1981_1990.nc"))
r_ps_91_00 <- rast(file.path(   
  pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_ps_global_daily_1991_2000.nc"))
r_ps_01_10 <- rast(file.path(   
  pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_ps_global_daily_2001_2010.nc"))
r_ps_11_14 <- rast(file.path(   
  pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_ps_global_daily_2011_2014.nc"))
r_ps <- c(r_ps_71_80, r_ps_81_90, r_ps_91_00, r_ps_01_10, r_ps_11_14)
time(r_ps) <- c(time(r_ps_71_80), time(r_ps_81_90), time(r_ps_91_00), time(r_ps_01_10), time(r_ps_11_14))

#longwave radiation (rlds)
r_rlds_71_80 <- rast(file.path(   
  pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_rlds_global_daily_1971_1980.nc"))
r_rlds_81_90 <- rast(file.path(   
  pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_rlds_global_daily_1981_1990.nc"))
r_rlds_91_00 <- rast(file.path(   
  pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_rlds_global_daily_1991_2000.nc"))
r_rlds_01_10 <- rast(file.path(   
  pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_rlds_global_daily_2001_2010.nc"))
r_rlds_11_14 <- rast(file.path(   
  pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_rlds_global_daily_2011_2014.nc"))
r_rlds <- c(r_rlds_71_80, r_rlds_81_90, r_rlds_91_00, r_rlds_01_10, r_rlds_11_14)
time(r_rlds) <- c(time(r_rlds_71_80), 
                  time(r_rlds_81_90), time(r_rlds_91_00), time(r_rlds_01_10), time(r_rlds_11_14))

#shortwave radiation (rsds)
r_rsds_71_80 <- rast(file.path(   
  pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_rsds_global_daily_1971_1980.nc"))
r_rsds_81_90 <- rast(file.path(   
  pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_rsds_global_daily_1981_1990.nc"))
r_rsds_91_00 <- rast(file.path(   
  pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_rsds_global_daily_1991_2000.nc"))
r_rsds_01_10 <- rast(file.path(   
  pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_rsds_global_daily_2001_2010.nc"))
r_rsds_11_14 <- rast(file.path(   
  pfad, "gfdl-esm4_r1i1p1f1_w5e5_historical_rsds_global_daily_2011_2014.nc"))
r_rsds <- c(r_rsds_71_80, r_rsds_81_90, r_rsds_91_00, r_rsds_01_10, r_rsds_11_14)
time(r_rsds) <- c(time(r_rsds_71_80), time(r_rsds_81_90), time(r_rsds_91_00), time(r_rsds_01_10), time(r_rsds_11_14))

#combine raster in one list 
raster <- list(
  r_tas  = r_tas,
  r_pr   = r_pr,
  r_ps   = r_ps,
  r_rlds = r_rlds,
  r_rsds = r_rsds
)

out_files <- lapply(
  X = EZG_id,
  FUN = raster_extraction,
  EZG_shape = EZG_shape,
  rasters   = raster,
  out_dir   = ordner_isimip_output, 
  scenario = scenario
)
