##02_02 Extraction of isimip raster data per catchment; creating of forcing data set for future period

#load packages 
library(terra)

#functions
source("R_func/raster_extraction.R")

#path to future ISIMIP raster data 
pfad <- "Data/ISIMIP_fut/"

#output folder
ordner_isimip_output <- "Output/ISIMIP/ISIMIP_forcing_future/ssp370"
scenario <- "ssp370"

#load catchment shapes and catchment ids 
EZG_shape <- vect("Data/Camels/camels-20250724T0232Z/basin_set_full_res/HCDN_nhru_final_671.shp")
EZG <- read.csv2("Data/Camels/camels-20250724T0232Z/camels_name.txt", header = TRUE, colClasses = c("character", "character", "character"))
EZG_id <- EZG[,1]

#load future raster data (ssp370)
# temperature (tas)
r_tas_41_50 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_tas_global_daily_2041_2050.nc"))
r_tas_51_60 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_tas_global_daily_2051_2060.nc"))
r_tas_61_70 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_tas_global_daily_2061_2070.nc"))
r_tas_71_80 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_tas_global_daily_2071_2080.nc"))
r_tas_81_90 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_tas_global_daily_2081_2090.nc"))
r_tas_91_00 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_tas_global_daily_2091_2100.nc"))
r_tas <- c(r_tas_41_50, r_tas_51_60, r_tas_61_70, r_tas_71_80, r_tas_81_90, r_tas_91_00)
time(r_tas) <- c(time(r_tas_41_50), time(r_tas_51_60), time(r_tas_61_70), time(r_tas_71_80), time(r_tas_81_90), time(r_tas_91_00))

# precipitation (pr)
r_pr_41_50 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_pr_global_daily_2041_2050.nc"))
r_pr_51_60 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_pr_global_daily_2051_2060.nc"))
r_pr_61_70 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_pr_global_daily_2061_2070.nc"))
r_pr_71_80 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_pr_global_daily_2071_2080.nc"))
r_pr_81_90 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_pr_global_daily_2081_2090.nc"))
r_pr_91_00 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_pr_global_daily_2091_2100.nc"))
r_pr <- c(r_pr_41_50, r_pr_51_60, r_pr_61_70, r_pr_71_80, r_pr_81_90, r_pr_91_00)
time(r_pr) <- c(time(r_pr_41_50), time(r_pr_51_60), time(r_pr_61_70), time(r_pr_71_80), time(r_pr_81_90), time(r_pr_91_00))

# air pressure (ps)
r_ps_41_50 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_ps_global_daily_2041_2050.nc"))
r_ps_51_60 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_ps_global_daily_2051_2060.nc"))
r_ps_61_70 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_ps_global_daily_2061_2070.nc"))
r_ps_71_80 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_ps_global_daily_2071_2080.nc"))
r_ps_81_90 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_ps_global_daily_2081_2090.nc"))
r_ps_91_00 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_ps_global_daily_2091_2100.nc"))
r_ps <- c(r_ps_41_50, r_ps_51_60, r_ps_61_70, r_ps_71_80, r_ps_81_90, r_ps_91_00)
time(r_ps) <- c(time(r_ps_41_50), time(r_ps_51_60), time(r_ps_61_70), time(r_ps_71_80), time(r_ps_81_90), time(r_ps_91_00))

# longwave radiation (rlds)
r_rlds_41_50 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_rlds_global_daily_2041_2050.nc"))
r_rlds_51_60 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_rlds_global_daily_2051_2060.nc"))
r_rlds_61_70 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_rlds_global_daily_2061_2070.nc"))
r_rlds_71_80 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_rlds_global_daily_2071_2080.nc"))
r_rlds_81_90 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_rlds_global_daily_2081_2090.nc"))
r_rlds_91_00 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_rlds_global_daily_2091_2100.nc"))
r_rlds <- c(r_rlds_41_50, r_rlds_51_60, r_rlds_61_70, r_rlds_71_80, r_rlds_81_90, r_rlds_91_00)
time(r_rlds) <- c(time(r_rlds_41_50), time(r_rlds_51_60), time(r_rlds_61_70), time(r_rlds_71_80), time(r_rlds_81_90), time(r_rlds_91_00))

# shortwave radiation (rsds)
r_rsds_41_50 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_rsds_global_daily_2041_2050.nc"))
r_rsds_51_60 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_rsds_global_daily_2051_2060.nc"))
r_rsds_61_70 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_rsds_global_daily_2061_2070.nc"))
r_rsds_71_80 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_rsds_global_daily_2071_2080.nc"))
r_rsds_81_90 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_rsds_global_daily_2081_2090.nc"))
r_rsds_91_00 <- rast(file.path(pfad, "gfdl-esm4_r1i1p1f1_w5e5_ssp370_rsds_global_daily_2091_2100.nc"))
r_rsds <- c(r_rsds_41_50, r_rsds_51_60, r_rsds_61_70, r_rsds_71_80, r_rsds_81_90, r_rsds_91_00)
time(r_rsds) <- c(time(r_rsds_41_50), time(r_rsds_51_60), time(r_rsds_61_70), time(r_rsds_71_80), time(r_rsds_81_90), time(r_rsds_91_00))

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

