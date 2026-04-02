calculate_changepoint <- function(topo_path, input_data, model_output,warmup, id, show_plot = TRUE){
  #compute soil moisture treshold based on changepoint method  
  #1) Identify independent runoff peak events (declustering) from Q_sim, using two conditions:
  # Q_sim > 10th percentile (potential event)
  # precipitation on the day before > 0
  #2) For each event, antecedent soil moisture (SM) is derived, at start of rainfall episode
  #3) If SM–Q relationship is significant (Spearman), exponential model is fitted and changepoint on the fitted curve detected (PELT)
  #input:
  #topo_path: path to camels_topo.txt 
  #input_data: model forcing including observed flow
  #model_output: hbv output containing qsim and states
  #warmup: number of warmup days to remove
  #id: catchment id 
  #show_plot: is plot generated?
  #output:
  #list with treshold_SM, runoff_peakd, nRMSE, has_cp (was a Changepoint detected), and plot
  
  p <- NULL #initialize plot object
  
  #extract HBV states from model output 
  storage_alle <- attr(model_output, "STATES")
  length_storage_df <- length(storage_alle[,1])
  n_dates <- nrow(input_data)
  
  #build base data frame and remove warm-up period
  base_df <- data.frame(date = input_data$date, Q_obs = input_data$flow, Q_sim = attr(model_output, "Q_sim"), prec = input_data$prec)
  base_df <- slice(base_df, (warmup + 1):nrow(base_df))
  
  #precipitation on the previous day
  base_df$prec_tag1 <- c(NA, base_df$prec[-nrow(base_df)])
  #define runoff treshold for potential event (10th percentile)
  runoff_10 <- quantile(base_df$Q_sim, 0.1, na.rm = TRUE)
  
  #define n_days
  camels_topo <- read.table(topo_path, header = TRUE, sep = ";")
  area <- camels_topo[camels_topo$gauge_id == as.numeric(id), 6]
  n_days <- round(5 + log(area))
  
  #mark potential peak (Q_sim above threshold and rainfall on previous day > 0)
  base_df$is_peak <- base_df$Q_sim > runoff_10 & base_df$prec_tag1 > 0
  
  #initialize event IDs for declustering
  base_df$event_id <- NA 
  peak_rows <- which(base_df$is_peak)
  
  event_id <- 1 
  group_indices <- c(peak_rows[1]) #indices of peaks belonging to the current event group

  #loop over all peaks (declustering)
  i <- 2
  while (i <= length(peak_rows)) {
    
    current_index <- peak_rows[i]
    
    # Peak with maximum Q within the current group (comparison peak Pi)
    group_q <- base_df$Q_sim[group_indices] 
    max_index <- group_indices[which.max(group_q)]
    
    #camparison with current peak
    #condition 1: time distance between peaks
    days_between <- as.numeric(base_df$date[current_index] - base_df$date[max_index])
    #condition 2: discharge between peaks must drop below 2/3 of the smaller peak
    between_rows <- which(base_df$date > base_df$date[max_index] &
                            base_df$date < base_df$date[current_index])
    
    min_q_between <- if (length(between_rows) == 0) Inf else min(base_df$Q_sim[between_rows], na.rm = TRUE)
    q_thresh <- (2/3) * min(base_df$Q_sim[current_index], base_df$Q_sim[max_index])
    
    if (days_between >= n_days && min_q_between < q_thresh) {
      # Both conditions satisfied -> finalize current event group;current event_id is assigned to the highest peak of the group
      base_df$event_id[max_index] <- event_id 
      event_id <- event_id + 1
      #start new group with current peak
      group_indices <- c(current_index) 
      i <- i + 1 
    } else {
      #conditions not fulfilled -> keep in same group
      group_indices <- c(group_indices, current_index)
      i <- i + 1
    }
  }
  
  #Finalize last group
  if (length(group_indices) > 0) {
    max_index <- group_indices[which.max(base_df$Q_sim[group_indices])]
    base_df$event_id[max_index] <- event_id
  }
  
  #extract independent runoff peak events 
  runoff_peaks <- base_df %>%
    filter(!is.na(event_id))
  
  #Determine start of rainfall episode before each runoff peak:
  #walk backwards from peak date until first day with 0 mm precipitation is found
  runoff_peaks$start_regenzeitraum <- as.Date(NA)
  
  for (i in 1:nrow(runoff_peaks)) {
    
    peak_date <- runoff_peaks$date[i]
    peak_pos <- which(base_df$date == peak_date)

    for (j in seq(peak_pos, 1, by = -1)) {
      if (base_df$prec[j] == 0) {
        runoff_peaks$start_regenzeitraum[i] <- base_df$date[j]
        break
      }
    }
  }
  
  #drop events where rainfall start could not be determined (outside record)
  runoff_peaks <- runoff_peaks %>%
    filter(!is.na(start_regenzeitraum))
  
  #extract soil moisture from HBV states
  if (length_storage_df == n_dates){
    Soil_moisture_df <- data.frame(date=input_data$date, soil_moisture = storage_alle[,3])
  } else {
    Soil_moisture_df <- data.frame(date=input_data$date, soil_moisture = storage_alle[2:length_storage_df,3])
  }
  
  #remove warump period
  Soil_moisture_df <- slice(Soil_moisture_df, (warmup + 1):nrow(Soil_moisture_df))
  
  #assign antecedent SM at the start of rainfall period 
  runoff_peaks$SM <- NA
  for (i in 1:nrow(runoff_peaks)) {
    SM_date <- runoff_peaks$start_regenzeitraum[i]
    SM_value <- Soil_moisture_df$soil_moisture[Soil_moisture_df$date == SM_date]
    runoff_peaks$SM[i] <- SM_value 
    
  }
  
  #built event data set with normalized SM values
  sm_min <- min(runoff_peaks$SM, na.rm = TRUE)
  sm_max <- max(runoff_peaks$SM, na.rm = TRUE)
  
  events_df <- data.frame(
    event_id = runoff_peaks$event_id,
    date = runoff_peaks$date, 
    Q_sim = runoff_peaks$Q_sim,
    SM_norm = (runoff_peaks$SM - sm_min) / (sm_max - sm_min)
  )
  
  #test spearman correlation (relationship runoff magnitude and antecedent soil moisture)
  spearman_cor <- cor.test(events_df$Q_sim, events_df$SM_norm, method = "spearman", exact = FALSE)
  p_value <- spearman_cor$p.value
  
  #if correlation significant, fit exponential SM-Q relationship and detect changpoint
  if (p_value < 0.05) {
    
    model <- try(
      nls(
        Q_sim ~ a * exp(b * SM_norm),
        data  = events_df,
        start = list(a = 1, b = 0.1)
      ),
      silent = TRUE
    )
    
  
    if (inherits(model, "try-error")){
      # model fit crasht
      print(paste0(id, "- model fit crasht"))
      return(list(
        threshold_SM = NA_real_,
        runoff_peaks = runoff_peaks,
        nRMSE        = NA_real_,
        has_cp       = FALSE,
        plot         = NULL
      ))
    }
    
    #changepoint detection on fitted Q(SM) curve (sorted by SM) (PELT)
    df_sorted <- events_df %>% arrange(SM_norm)
    df_sorted$Q_fit <- predict(model, newdata = data.frame(SM_norm = df_sorted$SM_norm))
    
    cpt <- cpt.mean(df_sorted$Q_fit, method = "PELT")
    change_index <- cpts(cpt)[1]
    threshold_SM <- df_sorted$SM_norm[change_index]
    
    #handle missing/invalid changepoint results
    if (is.na(threshold_SM)){ 
      print(paste0(id, ",Threshold = NA"))
      return(list(
        threshold_SM = NA_real_,
        runoff_peaks = runoff_peaks,
        nRMSE        = NA_real_,
        has_cp       = FALSE,
        plot         = NULL
      ))
    }
    
    #model fit metric: calculation of nRMSE
    y_obs <- df_sorted$Q_sim
    y_fit <- df_sorted$Q_fit
    res <- y_obs - y_fit
    nRMSE <- sqrt(mean(res^2)) / sd(y_obs) 
    
    #plot 
    if (show_plot) {
      fit_curve <- data.frame(
        SM_norm = seq(min(df_sorted$SM_norm, na.rm = TRUE),
                      max(df_sorted$SM_norm, na.rm = TRUE),
                      length.out = 200)
      )
      fit_curve$Q_fit <- predict(model, newdata = fit_curve)
      
      # p-Text (mit "< 0.001")
      p_val <- as.numeric(spearman_cor$p.value)
      p_txt <- if (is.finite(p_val) && p_val < 0.001) "< 0.001" else sprintf("%.3f", p_val)
      
      y_max <- max(df_sorted$Q_sim, na.rm = TRUE)
      
      p <- ggplot(df_sorted, aes(x = SM_norm, y = Q_sim)) +
        geom_point(shape = 21, fill = NA, color = "blue", stroke = 0.8, size = 2.2) +
        geom_line(data = fit_curve, aes(SM_norm, Q_fit), color = "darkorange",
                  linetype = "dashed", linewidth = 1) +
        geom_vline(xintercept = threshold_SM, color = "red", linetype = "dashed") +
        annotate("text", x = 0.1, y = 0.78 * y_max,
                 label = paste0("ChangePoint = ", round(threshold_SM, 2)),
                 hjust = 0) +
        annotate("text", x = 0.1, y = 0.70 * y_max,
                 label = paste0("p = ", p_txt),
                 hjust = 0) +
        labs(
          title = paste0("SM–Q Relationship with Change Point: ", id),
          x = "Soil moisture [0–1]",
          y = expression("Event runoff [m"^3*"/s]")
        ) +
        theme_minimal() +
        expand_limits(x = 0, y = 0)
      
      # #plot: example catchement
      # p2 <- ggplot(df_sorted, aes(x = SM_norm, y = Q_sim)) +
      #   geom_point(shape = 21, fill = NA, color = "blue", stroke = 0.8, size = 2.2, alpha = 0.4) +
      #   geom_line(data = fit_curve, aes(SM_norm, Q_fit), color = "darkorange",
      #             linetype = "dashed", linewidth = 0.8) +
      #   geom_vline(xintercept = threshold_SM, color = "red", linetype = "dashed",linewidth = 0.6, alpha = 0.8 ) } +
      #   annotate("text", x = 0.1, y = 0.78 * y_max,
      #            label = paste0("Changepoint = ", round(threshold_SM, 2)),
      #            hjust = 0) +
      #   # annotate("text", x = 0.1, y = 0.70 * y_max,
      #   #          label = paste0("p = ", p_txt),
      #   #          hjust = 0) +
      #   labs(
      #     title = paste0("Soil moisture–runoff relationship with identified changepoint (example catchment)"),
      #     x = "Antecedent soil moisture [0–1]",
      #     y = expression("Event runoff [m"^3*"/s]")
      #   ) +
      #   theme_minimal() +
      #   expand_limits(x = 0, y = 0)
      # 
      
      return(list(threshold_SM = threshold_SM, runoff_peaks = runoff_peaks, nRMSE = nRMSE, has_cp = TRUE, plot = p))
    } 
    
    # show_plot = FALSE
    return(list(threshold_SM = threshold_SM,
                runoff_peaks = runoff_peaks,
                nRMSE = nRMSE,
                has_cp = TRUE,
                plot = NULL))
  } else {
    
    #no significant correlation  
    print(paste0(id, " - keine signifikante Korrelation"))
    return(list(
      threshold_SM = NA_real_,
      runoff_peaks = runoff_peaks,
      nRMSE        = NA_real_,
      has_cp       = FALSE,
      plot         = NULL
    ))
    
  }
} 

