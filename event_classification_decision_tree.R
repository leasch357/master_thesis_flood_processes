##Funktion für Decision tree 
event_classification_quantile = function(decision_df, parts_rainsnow = 1/3, parts_fracextreme = 2/3, sat_thresh = 0.9){
  #classification tree for flood peak classification
  #Input:
  #decision df: decision data frame generated from event_classification_df function
  #parts_rainsnow: threshold limit overlap between rainfall and snowmelt
  #parts_fracextreme: threshold limit for how much weekly rainfall has to fall on one day to be considered a short rainfall
  #sat_thresh: threshold limit saturated soil
  #Output: 
  #decision df with additional flood process column
  
  ELSE = TRUE
  decision_result = decision_df %>% mutate(.,mech_out = with(.,case_when(
    #minimum rainfall limit
    (Ptotal < 1) ~ "noclass",
    (!is.na(S7) & (S7/Ptotal)>=parts_rainsnow & (P7/Ptotal)>=parts_rainsnow) ~ "rainandsnow",
    (!is.na(S7) & S7>=S7_quant) ~ "snowmelt",
    (sat_start >= sat_thresh & P7 >= P7_mean) ~ "soilsat",
    (P7 >= P7_quant & Pmax_frac >= parts_fracextreme) ~ "rainfall",
    (P7 >= P7_quant & Pmax_frac < parts_fracextreme) ~ "longrainfall",
    (Pmax >= Pday_quant) ~ "rainfall",
    (Pmax < Pday_quant) ~ "noclass",
    ELSE ~ "missingData"
  )))
  
  #turn numbers back into dates
  decision_result = decision_result %>% mutate(flood_date = as.Date(flood_date,  origin = "1970-01-01"))
  return(decision_result)
}
