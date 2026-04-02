amax <- function(input_data, model_output,warmup){
  #calculate amax events
  #input
  #input_data Modelinput data
  #model_output 
  #warmup
  #output:
  #data frame with flood date and magnnitude 

  output_data <- data.frame(date = input_data$date, Q_sim = attr(model_output, "Q_sim"))
  output_data <- slice(output_data, (warmup + 1):nrow(output_data))
  
  flood_data_df <- output_data %>%
    mutate(jahr= year(date)) %>%
    group_by(jahr) %>% 
    slice_max(order_by = Q_sim, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  flood_data_df <- flood_data_df[, c("date", "Q_sim")]
  colnames(flood_data_df) = c("flood_date", "flood_magn")
  
  return(flood_data_df)
}

