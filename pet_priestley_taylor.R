pet_priestley_taylor <- function(forcing_df, alpha = 1.26, albedo  = 0.23, emis    = 0.95){
  #compute daily PET values 
  #input: 
  #forcing_df: data.frame with date, temp (C°), prec (mm), ps (Kpa), rlds (W m-2), rsds (W m-2)
  #alpha = 1.26, standart alpha value; albedo  = 0.23, standart albedo following FAO; emis    = 0.95 standart emissivity following FAO
  #output: 
  #pet df with daily PET values 
  
  #constants
  #G <- 0 soil heat flux 
  sigma <- 5.670374e-8      # W m-2 K-4
  cp <- 1.013e-3         # MJ kg-1 °C-1
  eps   <- 0.622            # Mw ratio of water vapor/dry air
  Wm2_to_MJ <- 86400 / 1e6  #1 W/m2 = 1 J/m2 s
  
  pet_df <- forcing_df
  
  #FAO 
  pet_df$lambda = 2.501 - (0.002361 * pet_df$temp) # MJ kg-1
  pet_df$gamma = (cp * pet_df$ps) / (eps * pet_df$lambda) # kPa °C-1
  pet_df$esat <- 0.6108 * exp(17.27 * pet_df$temp / (pet_df$temp + 237.3)) # kPa
  pet_df$s <- 4098 * pet_df$esat / (pet_df$temp + 237.3)^2 # kPa °C-1
  
  pet_df$Rns <- (1-albedo) * pet_df$rsds  
  pet_df$Rnl <- pet_df$rlds - emis *sigma* ((pet_df$temp + 273.15)^4) 
  pet_df$Rn <- (pet_df$Rns + pet_df$Rnl) * Wm2_to_MJ
  
  pet_df$PET <- alpha/pet_df$lambda * pet_df$s/(pet_df$s+pet_df$gamma) * pet_df$Rn #PET in mm/d 
  pet_df$PET[pet_df$PET < 0] <- 0
  
  return(pet_df = pet_df)
}