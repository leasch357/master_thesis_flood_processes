get_rel_months <- function(df, date_col = "date") {
  #calculate monthly frequencies of the amax events
  #input:
  #df dataframe with magnitude and date; datecol  
  dates <- as.Date(df[[date_col]])
  months <- as.integer(format(dates, "%m"))
  rel_months <- prop.table(table(factor(months, levels = 1:12)))
  as.numeric(rel_months)
}

doy_to_label <- function(doy) {
  #transfers doy back into date 
  ref_date <- as.Date("2021-01-01") + round(doy) - 1
  format(ref_date, "%d.%m")
}

circ_stats_doy <- function(date_col, d = 365) {
  #circular statistics to calculate mean occurence day, using circ_stats 
  #input:
  #date_col column with date of amax events 
  #output:
  #list with mean day of occurence (doy, day) and R
  
  # transfer in doy
  doy_vec <- yday(date_col)
  
  #convert to angles
  ang <- doy_vec * 2 * pi / d 
  
  #calculate mean x and y coordinate
  x_bar <- mean(cos(ang)) 
  y_bar <- mean(sin(ang)) 
  
  #concentration R from vector length 
  R     <- sqrt(x_bar^2 + y_bar^2) 
  
  # mean angle
  alpha <- atan2(y_bar, x_bar) 
  if (alpha < 0) alpha <- alpha + 2 * pi
  
  #convert back to doy
  mean_doy <- alpha * d / (2 * pi) 
  
  list(
    mean_doy   = mean_doy,
    mean_day = doy_to_label(mean_doy),
    R          = R
  )
}

circ_diff_days <- function(d1, d2, d = 365) {
  #calculate the total difference between two days
  dif <- round(abs(d1 - d2))
  pmin(dif, d - dif)
}

