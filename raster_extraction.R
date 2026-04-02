raster_extraction <- function(id, EZG_shape, rasters, out_dir, scenario){
  #extract catchment mean isimip forcing time series and write forcing csv
  #input: 
  #id: catchment id 
  #ezg shape: catchment polygons (terra vect)
  #rasters: list with isimip raster stacks (r tas,r pr,r ps,r rlds,r rsds)
  #out dir: output directory for forcing csv files
  #scenario: scenario label used in output file name (e.g. hist,ssp370,ssp585)
  #output:
  #returns path to written forcing csv
  
  #select catchment and convert catchments CRS to grid CRS
  ezg <- EZG_shape[EZG_shape$hru_id == as.numeric(id), ]
  ezg <- project(ezg, crs(rasters$r_tas))

  #extract dates from r_tas 
  dates <- as.Date(time(rasters$r_tas))
  
  #extract are mean 
  tas_ext  <- extract(rasters$r_tas,  ezg, fun = mean, na.rm = TRUE, exact = TRUE)
  pr_ext   <- extract(rasters$r_pr,   ezg, fun = mean, na.rm = TRUE, exact = TRUE)
  ps_ext   <- extract(rasters$r_ps,   ezg, fun = mean, na.rm = TRUE, exact = TRUE)
  rlds_ext <- extract(rasters$r_rlds, ezg, fun = mean, na.rm = TRUE, exact = TRUE)
  rsds_ext <- extract(rasters$r_rsds, ezg, fun = mean, na.rm = TRUE, exact = TRUE)
  
  tas_K <- as.numeric(tas_ext[1, -1])
  pr    <- as.numeric(pr_ext[1,  -1])
  ps    <- as.numeric(ps_ext[1,  -1])
  rlds  <- as.numeric(rlds_ext[1, -1])
  rsds  <- as.numeric(rsds_ext[1, -1])
  
  #convert units and build forcing data frame
  forcing_df <- data.frame(
    date = dates,
    temp = tas_K - 273.15,  # K -> °C
    prec = pr * 86400,      # kg m-2 s-1 -> mm d-1
    ps   = ps / 1000,       # Pa -> kPa
    rlds = rlds,            # W m-2 
    rsds = rsds             # W m-2
  )
  
  csv_path <- file.path(out_dir, paste0("Forcing_df_", id,"_", scenario, ".csv"))
  write.csv(forcing_df, csv_path, row.names = FALSE)
  
  return(csv_path)
}
