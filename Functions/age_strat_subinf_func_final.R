# version 4. DALY model 

age_subinf_adult_psa <- function(infection_df, lhs_sample) {
  
  result_list <- vector("list", length = 16)  # Initialize a list to store results
  
  df_symp <- df_asymp <- df_hosp <- df_hosp_acute <- df_hosp_chronic_6m <- df_hosp_chronic_12m <- df_hosp_chronic_30m <- df_hosp_subac <- df_hosp_fatal <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
  colnames(df_symp) <- colnames(df_hosp) <- colnames(df_asymp) <- colnames(df_hosp_acute) <- colnames(df_hosp_chronic_6m) <- colnames(df_hosp_chronic_12m) <- colnames(df_hosp_chronic_30m) <- colnames(df_hosp_subac) <- colnames(df_hosp_fatal) <- paste0("Run_", 1:nrow(lhs_sample))
  
  df_nonhosp <- df_nonhosp_acute <- df_nonhosp_chronic_6m <- df_nonhosp_chronic_12m <- df_nonhosp_chronic_30m <- df_nonhosp_subac <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
  colnames(df_nonhosp) <- colnames(df_nonhosp_acute) <- colnames(df_nonhosp_chronic_6m) <- colnames(df_nonhosp_chronic_12m) <- colnames(df_nonhosp_chronic_30m) <- colnames(df_nonhosp_subac) <- paste0("Run_", 1:nrow(lhs_sample))
  
  df_yld_hosp_chronic <- df_yld_hosp_acute <- df_yld_hosp_subac <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
  colnames(df_yld_hosp_chronic) <- colnames(df_yld_hosp_acute) <- colnames(df_yld_hosp_subac) <- paste0("Run_", 1:nrow(lhs_sample))
  
  df_yld_nh_chronic <- df_yld_nh_acute <- df_yld_nh_subac <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
  colnames(df_yld_nh_chronic) <- colnames(df_yld_nh_acute) <- colnames(df_yld_nh_subac) <- paste0("Run_", 1:nrow(lhs_sample))
  
  df_yld_acute <- df_yld_subac <- df_yld_chronic <- df_yll <- df_daly <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
  colnames(df_yld_acute) <- colnames(df_yld_subac) <- colnames(df_yld_chronic) <- colnames(df_yll) <- colnames(df_daly) <- paste0("Run_", 1:nrow(lhs_sample))
  
  
  for(k in 1:nrow(lhs_sample)) {
    
    for(i in 1:nrow(infection_df)) {
      
      # Mid estimates
      med_inf = infection_df$tot_infec_med_eld[i]
      
      # Apply Asian symptomatic to Asian continent 
      
      if (infection_df$continent[i] == "Asia") {
        
        df_symp[i, k] <- med_inf * lhs_sample$symp_asia[k]
        
        # Apply ECSA symptomatic to African continent
        
      } else if (infection_df$continent[i] == "Africa") {
        
        df_symp[i, k]  <- med_inf * lhs_sample$symp_africa[k]
        
        # Apply ECSA symptomatic to African continent
        
      } else if (infection_df$continent[i] == "North America") {
        
        df_symp[i, k]  <- med_inf * lhs_sample$symp_america[k]
        
      } else if (infection_df$continent[i] == "South America") {
        
        df_symp[i, k]  <- med_inf * lhs_sample$symp_america[k]
        
        # Apply overall symptomatic to other continents (EU/Oceania)
        
      } else {
        
        df_symp[i, k] <- med_inf * lhs_sample$symp_overall[k]
        
      }
      
      
      df_asymp[i, k]               <- med_inf - df_symp[i, k]
      df_hosp[i, k]                <- df_symp[i, k]*lhs_sample$hosp[k]
      df_hosp_acute[i, k]          <- df_hosp[i, k]*lhs_sample$acute[k]
      df_hosp_chronic_6m[i, k]     <- df_hosp[i, k]*lhs_sample$chr6m[k] 
      df_hosp_chronic_12m[i, k]    <- df_hosp[i, k]*lhs_sample$chr12m[k] 
      df_hosp_chronic_30m[i, k]    <- df_hosp[i, k]*lhs_sample$chr30m[k]
      df_hosp_subac[i, k]          <- df_hosp[i, k]*lhs_sample$subac[k]

      df_nonhosp[i, k]              <- df_symp[i, k] - df_hosp[i, k]
      df_nonhosp_acute[i, k]        <- df_nonhosp[i, k]*lhs_sample$acute[k]
      df_nonhosp_chronic_6m[i ,k]   <- df_nonhosp[i, k]*lhs_sample$chr6m[k]
      df_nonhosp_chronic_12m[i ,k]  <- df_nonhosp[i, k]*lhs_sample$chr12m[k] 
      df_nonhosp_chronic_30m[i ,k]  <- df_nonhosp[i, k]*lhs_sample$chr30m[k]
      df_nonhosp_subac[i, k]        <- df_nonhosp[i, k]*lhs_sample$subac[k]
      df_hosp_fatal[i, k]           <- df_hosp[i, k]*lhs_sample$fatal_hosp[k] + df_nonhosp[i, k]*lhs_sample$fatal_nonhosp[k]

      # YLD estimation
      
      df_yld_hosp_chronic[i, k]        <- df_hosp_chronic_6m[i, k]*lhs_sample$dur_6m[k]*lhs_sample$dw_chronic[k] + df_hosp_chronic_12m[i, k]*lhs_sample$dur_12m[k]*lhs_sample$dw_chronic[k] + df_hosp_chronic_30m[i, k]*lhs_sample$dur_30m[k]*lhs_sample$dw_chronic[k]
      df_yld_hosp_acute[i, k]          <- df_hosp_acute[i, k]*lhs_sample$dw_hosp[k]*lhs_sample$dur_acute[k]
      df_yld_hosp_subac[i, k]          <- df_hosp_subac[i, k]*lhs_sample$dw_subac[k]*lhs_sample$dur_subac[k]

      df_yld_nh_chronic[i, k]      <- df_nonhosp_chronic_6m[i, k]*lhs_sample$dur_6m[k]*lhs_sample$dw_chronic[k] + df_nonhosp_chronic_12m[i, k]*lhs_sample$dur_12m[k]*lhs_sample$dw_chronic[k] + df_nonhosp_chronic_30m[i, k]*lhs_sample$dur_30m[k]*lhs_sample$dw_chronic[k]
      df_yld_nh_acute[i, k]        <- df_nonhosp_acute[i, k]*lhs_sample$dw_hosp[k]*lhs_sample$dur_acute[k]
      df_yld_nh_subac[i, k]        <- df_nonhosp_subac[i, k]*lhs_sample$dw_subac[k]*lhs_sample$dur_subac[k]

      # final daly 
      
      df_yld_acute[i, k]                <- df_yld_hosp_acute[i, k] + df_yld_nh_acute[i, k] 
      df_yld_subac[i, k]                <- df_yld_hosp_subac[i, k] + df_yld_nh_subac[i, k] 
      df_yld_chronic[i, k]              <- df_yld_hosp_chronic[i, k] + df_yld_nh_chronic[i, k] 
      df_yll[i, k]                      <- df_hosp_fatal[i, k] * lhs_sample$le_lost[k]
      df_daly[i, k]                     <- df_yld_acute[i, k] + df_yld_subac[i, k] + df_yld_chronic[i, k] + df_yll[i, k] 
        
    }
  }
  
  result_list <- list(symp                   = df_symp, 
                      asymp                  = df_asymp, 
                      df_hosp                = df_hosp, 
                      df_hosp_acute          = df_hosp_acute,
                      df_hosp_chronic_6m     = df_hosp_chronic_6m,
                      df_hosp_chronic_12m    = df_hosp_chronic_12m,
                      df_hosp_chronic_30m    = df_hosp_chronic_30m, 
                      df_hosp_subac          = df_hosp_subac, 
                      df_nonhosp             = df_nonhosp, 
                      df_nonhosp_acute       = df_nonhosp_acute,
                      df_nonhosp_chronic_6m  = df_nonhosp_chronic_6m, 
                      df_nonhosp_chronic_12m = df_nonhosp_chronic_12m,
                      df_nonhosp_chronic_30m = df_nonhosp_chronic_30m,
                      df_nonhosp_subac       = df_nonhosp_subac,
                      df_yld_hosp_chronic    = df_yld_hosp_chronic, 
                      df_yld_hosp_acute      = df_yld_hosp_acute, 
                      df_yld_hosp_subac      = df_yld_hosp_subac,
                      df_yld_nh_chronic      = df_yld_nh_chronic, 
                      df_yld_nh_acute        = df_yld_nh_acute, 
                      df_yld_nh_subac        = df_yld_nh_subac,
                      df_yld_acute           = df_yld_acute, 
                      df_yld_subac           = df_yld_subac, 
                      df_yld_chronic         = df_yld_chronic, 
                      df_yll                 = df_yll, 
                      df_daly                = df_daly)  
  
  for (j in names(result_list)) {
    result_list[[j]] <- cbind(infection_df[, c( "country", "continent")], result_list[[j]])
    result_list[[j]]$mid <- apply(result_list[[j]][, c(3:1002)], 1, median, na.rm = TRUE)
    result_list[[j]]$lo  <- apply(result_list[[j]][, c(3:1002)], 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
    result_list[[j]]$hi  <- apply(result_list[[j]][, c(3:1002)], 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))
  }
  return(result_list)
}

age_subinf_child_psa <- function(infection_df, lhs_sample) {
  
  result_list <- vector("list", length = 16)  # Initialize a list to store results
  
  df_symp <- df_asymp <- df_hosp <- df_hosp_acute <- df_hosp_chronic_6m <- df_hosp_chronic_12m <- df_hosp_chronic_30m <- df_hosp_subac <- df_hosp_fatal <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
  colnames(df_symp) <- colnames(df_hosp) <- colnames(df_asymp) <- colnames(df_hosp_acute) <- colnames(df_hosp_chronic_6m) <- colnames(df_hosp_chronic_12m) <- colnames(df_hosp_chronic_30m) <- colnames(df_hosp_subac) <- colnames(df_hosp_fatal) <- paste0("Run_", 1:nrow(lhs_sample))
  
  df_nonhosp <- df_nonhosp_acute <- df_nonhosp_chronic_6m <- df_nonhosp_chronic_12m <- df_nonhosp_chronic_30m <- df_nonhosp_subac <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
  colnames(df_nonhosp) <- colnames(df_nonhosp_acute) <- colnames(df_nonhosp_chronic_6m) <- colnames(df_nonhosp_chronic_12m) <- colnames(df_nonhosp_chronic_30m) <- colnames(df_nonhosp_subac) <- paste0("Run_", 1:nrow(lhs_sample))
  
  df_yld_hosp_chronic <- df_yld_hosp_acute <- df_yld_hosp_subac <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
  colnames(df_yld_hosp_chronic) <- colnames(df_yld_hosp_acute) <- colnames(df_yld_hosp_subac) <- paste0("Run_", 1:nrow(lhs_sample))
  
  df_yld_nh_chronic <- df_yld_nh_acute <- df_yld_nh_subac <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
  colnames(df_yld_nh_chronic) <- colnames(df_yld_nh_acute) <- colnames(df_yld_nh_subac) <- paste0("Run_", 1:nrow(lhs_sample))
  
  df_yld_acute <- df_yld_subac <- df_yld_chronic <- df_yll <- df_daly <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
  colnames(df_yld_acute) <- colnames(df_yld_subac) <- colnames(df_yld_chronic) <- colnames(df_yll) <- colnames(df_daly) <- paste0("Run_", 1:nrow(lhs_sample))
  
  
  for(k in 1:nrow(lhs_sample)) {
    
    for(i in 1:nrow(infection_df)) {
      
      # Mid estimates
      med_inf = infection_df$tot_infec_med_child[i]
      
      # Apply Asian symptomatic to Asian continent 
      
      if (infection_df$continent[i] == "Asia") {
        
        df_symp[i, k] <- med_inf * lhs_sample$symp_asia[k]
        
        # Apply ECSA symptomatic to African continent
        
      } else if (infection_df$continent[i] == "Africa") {
        
        df_symp[i, k]  <- med_inf * lhs_sample$symp_africa[k]
        
        # Apply ECSA symptomatic to African continent
        
      } else if (infection_df$continent[i] == "North America") {
        
        df_symp[i, k]  <- med_inf * lhs_sample$symp_america[k]
        
      } else if (infection_df$continent[i] == "South America") {
        
        df_symp[i, k]  <- med_inf * lhs_sample$symp_america[k]
        
        # Apply overall symptomatic to other continents (EU/Oceania)
        
      } else {
        
        df_symp[i, k] <- med_inf * lhs_sample$symp_overall[k]
        
      }
      
      
      df_asymp[i, k]               <- med_inf - df_symp[i, k]
      df_hosp[i, k]                <- df_symp[i, k]*lhs_sample$hosp[k]
      df_hosp_acute[i, k]          <- df_hosp[i, k]*lhs_sample$acute[k]
      df_hosp_chronic_6m[i, k]     <- df_hosp[i, k]*lhs_sample$chr6m[k] 
      df_hosp_chronic_12m[i, k]    <- df_hosp[i, k]*lhs_sample$chr12m[k] 
      df_hosp_chronic_30m[i, k]    <- df_hosp[i, k]*lhs_sample$chr30m[k]
      df_hosp_subac[i, k]          <- df_hosp[i, k]*lhs_sample$subac[k]

      df_nonhosp[i, k]              <- df_symp[i, k] - df_hosp[i, k]
      df_nonhosp_acute[i, k]        <- df_nonhosp[i, k]*lhs_sample$acute[k]
      df_nonhosp_chronic_6m[i ,k]   <- df_nonhosp[i, k]*lhs_sample$chr6m[k]
      df_nonhosp_chronic_12m[i ,k]  <- df_nonhosp[i, k]*lhs_sample$chr12m[k] 
      df_nonhosp_chronic_30m[i ,k]  <- df_nonhosp[i, k]*lhs_sample$chr30m[k]
      df_nonhosp_subac[i, k]        <- df_nonhosp[i, k]*lhs_sample$subac[k]
      df_hosp_fatal[i, k]           <- df_hosp[i, k]*lhs_sample$fatal_hosp[k] + df_nonhosp[i, k]*lhs_sample$fatal_nonhosp[k]
      
      
      # YLD estimation
      
      df_yld_hosp_chronic[i, k]        <- df_hosp_chronic_6m[i, k]*lhs_sample$dur_6m[k]*lhs_sample$dw_chronic[k] + df_hosp_chronic_12m[i, k]*lhs_sample$dur_12m[k]*lhs_sample$dw_chronic[k] + df_hosp_chronic_30m[i, k]*lhs_sample$dur_30m[k]*lhs_sample$dw_chronic[k]
      df_yld_hosp_acute[i, k]          <- df_hosp_acute[i, k]*lhs_sample$dw_hosp[k]*lhs_sample$dur_acute[k]
      df_yld_hosp_subac[i, k]          <- df_hosp_subac[i, k]*lhs_sample$dw_subac[k]*lhs_sample$dur_subac[k]
      
      df_yld_nh_chronic[i, k]      <- df_nonhosp_chronic_6m[i, k]*lhs_sample$dur_6m[k]*lhs_sample$dw_chronic[k] + df_nonhosp_chronic_12m[i, k]*lhs_sample$dur_12m[k]*lhs_sample$dw_chronic[k] + df_nonhosp_chronic_30m[i, k]*lhs_sample$dur_30m[k]*lhs_sample$dw_chronic[k]
      df_yld_nh_acute[i, k]        <- df_nonhosp_acute[i, k]*lhs_sample$dw_hosp[k]*lhs_sample$dur_acute[k]
      df_yld_nh_subac[i, k]        <- df_nonhosp_subac[i, k]*lhs_sample$dw_subac[k]*lhs_sample$dur_subac[k]
      
      # final daly 
      
      df_yld_acute[i, k]                <- df_yld_hosp_acute[i, k] + df_yld_nh_acute[i, k] 
      df_yld_subac[i, k]                <- df_yld_hosp_subac[i, k] + df_yld_nh_subac[i, k] 
      df_yld_chronic[i, k]              <- df_yld_hosp_chronic[i, k] + df_yld_nh_chronic[i, k] 
      df_yll[i, k]                      <- df_hosp_fatal[i, k] * lhs_sample$le_lost[k]
      df_daly[i, k]                     <- df_yld_acute[i, k] + df_yld_subac[i, k] + df_yld_chronic[i, k] + df_yll[i, k] 
      
    }
  }
  
  result_list <- list(symp                   = df_symp, 
                      asymp                  = df_asymp, 
                      df_hosp                = df_hosp, 
                      df_hosp_acute          = df_hosp_acute,
                      df_hosp_chronic_6m     = df_hosp_chronic_6m,
                      df_hosp_chronic_12m    = df_hosp_chronic_12m,
                      df_hosp_chronic_30m    = df_hosp_chronic_30m, 
                      df_hosp_subac          = df_hosp_subac, 
                      df_nonhosp             = df_nonhosp, 
                      df_nonhosp_acute       = df_nonhosp_acute,
                      df_nonhosp_chronic_6m  = df_nonhosp_chronic_6m, 
                      df_nonhosp_chronic_12m = df_nonhosp_chronic_12m,
                      df_nonhosp_chronic_30m = df_nonhosp_chronic_30m,
                      df_nonhosp_subac       = df_nonhosp_subac,
                      df_yld_hosp_chronic    = df_yld_hosp_chronic, 
                      df_yld_hosp_acute      = df_yld_hosp_acute, 
                      df_yld_hosp_subac      = df_yld_hosp_subac,
                      df_yld_nh_chronic      = df_yld_nh_chronic, 
                      df_yld_nh_acute        = df_yld_nh_acute, 
                      df_yld_nh_subac        = df_yld_nh_subac,
                      df_yld_acute           = df_yld_acute, 
                      df_yld_subac           = df_yld_subac, 
                      df_yld_chronic         = df_yld_chronic, 
                      df_yll                 = df_yll, 
                      df_daly                = df_daly)  
  
  for (j in names(result_list)) {
    result_list[[j]] <- cbind(infection_df[, c( "country", "continent")], result_list[[j]])
    result_list[[j]]$mid <- apply(result_list[[j]][, c(3:1002)], 1, median, na.rm = TRUE)
    result_list[[j]]$lo  <- apply(result_list[[j]][, c(3:1002)], 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
    result_list[[j]]$hi  <- apply(result_list[[j]][, c(3:1002)], 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))
  }
  return(result_list)
}


age_subinf_child_psa <- function(infection_df, lhs_sample) {
  
  result_list <- vector("list", length = 16)  # Initialize a list to store results
  
  df_symp <- df_asymp <- df_hosp <- df_hosp_recov <- df_hosp_chronic  <- df_hosp_chronic_mild <- df_hosp_chronic_severe <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
  colnames(df_symp) <- colnames(df_asymp) <- colnames(df_hosp) <- colnames(df_hosp_recov) <- colnames(df_hosp_chronic) <- colnames(df_hosp_chronic_mild) <- colnames(df_hosp_chronic_severe) <- paste0("Run_", 1:nrow(lhs_sample))
  
  df_nonhosp_severe <- df_nonhosp_severe_recov <- df_nonhosp_severe_chronic <- df_nonhosp_severe_chronic_mild <- df_nonhosp_severe_chronic_severe <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
  colnames(df_nonhosp_severe) <- colnames(df_nonhosp_severe_recov) <- colnames(df_nonhosp_severe_chronic)  <- colnames(df_nonhosp_severe_chronic_mild) <- colnames(df_nonhosp_severe_chronic_severe) <- paste0("Run_", 1:nrow(lhs_sample))
  
  df_nonhosp_mod <- df_nonhosp_mod_recov <- df_nonhosp_mod_chronic  <- df_nonhosp_mod_chronic_mild <- df_nonhosp_mod_chronic_severe <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
  colnames(df_nonhosp_mod) <- colnames(df_nonhosp_mod_recov) <- colnames(df_nonhosp_mod_chronic)  <- colnames(df_nonhosp_mod_chronic_mild) <- colnames(df_nonhosp_mod_chronic_severe) <- paste0("Run_", 1:nrow(lhs_sample))
  
  df_nonhosp_mild <- df_nonhosp_mild_recov <- df_nonhosp_mild_chronic  <- df_nonhosp_mild_chronic_mild <- df_hosp_fatal <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
  colnames(df_nonhosp_mild) <- colnames(df_nonhosp_mild_recov) <- colnames(df_nonhosp_mild_chronic)  <- colnames(df_nonhosp_mild_chronic_mild) <- colnames(df_hosp_fatal) <- paste0("Run_", 1:nrow(lhs_sample))
  
  df_yld_hosp_chronic <- df_yld_hosp_acute  <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
  colnames(df_yld_hosp_chronic) <- colnames(df_yld_hosp_acute)  <- paste0("Run_", 1:nrow(lhs_sample))
  
  df_yld_nh_severe_chronic <- df_yld_nh_severe_acute  <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
  colnames(df_yld_nh_severe_chronic) <- colnames(df_yld_nh_severe_acute)  <- paste0("Run_", 1:nrow(lhs_sample))
  
  df_yld_nh_mod_chronic <- df_yld_nh_mod_acute <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
  colnames(df_yld_nh_mod_chronic) <- colnames(df_yld_nh_mod_acute) <- paste0("Run_", 1:nrow(lhs_sample))
  
  df_yld_nh_mild_chronic <- df_yld_nh_mild_acute  <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
  colnames(df_yld_nh_mild_chronic) <- colnames(df_yld_nh_mild_acute)  <- paste0("Run_", 1:nrow(lhs_sample))
  
  df_yld_acute  <- df_yld_chronic <- df_yll <- df_daly <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
  colnames(df_yld_acute)  <- colnames(df_yld_chronic) <- colnames(df_yll) <- colnames(df_daly) <- paste0("Run_", 1:nrow(lhs_sample))
  
  
  for(k in 1:nrow(lhs_sample)) {
    
    for(i in 1:nrow(infection_df)) {
      
      # Mid estimates
      med_inf = infection_df$tot_infec_med_child[i]
      
      # Apply Asian symptomatic to Asian continent 
      
      if (infection_df$continent[i] == "Asia") {
        
        df_symp[i, k] <- med_inf * lhs_sample$symp_asia[k]
        
        # Apply ECSA symptomatic to African continent
        
      } else if (infection_df$continent[i] == "Africa") {
        
        df_symp[i, k]  <- med_inf * lhs_sample$symp_africa[k]
        
        # Apply ECSA symptomatic to African continent
        
      } else if (infection_df$continent[i] == "North America") {
        
        df_symp[i, k]  <- med_inf * lhs_sample$symp_america[k]
        
      } else if (infection_df$continent[i] == "South America") {
        
        df_symp[i, k]  <- med_inf * lhs_sample$symp_america[k]
        
        # Apply overall symptomatic to other continents (EU/Oceania)
        
      } else {
        
        df_symp[i, k] <- med_inf * lhs_sample$symp_overall[k]
        
      }
      
      
      df_asymp[i, k]               <- med_inf - df_symp[i, k]
      df_hosp[i, k]                <- df_symp[i, k]*lhs_sample$hosp[k]
      df_hosp_chronic[i, k]        <- df_hosp[i, k]*lhs_sample$severe_chronic[k]
      df_hosp_recov[i, k]          <- df_hosp[i, k]*(1 - (lhs_sample$severe_chronic[k]))
      df_hosp_chronic_mild[i, k]   <- df_hosp_chronic[i, k]*lhs_sample$severe_mild[k]
      df_hosp_chronic_severe[i, k] <- df_hosp_chronic[i, k]*lhs_sample$severe_severe[k]
      
      df_nonhosp_severe[i, k]                <- df_symp[i, k]*lhs_sample$nonhosp_severe[k]
      df_nonhosp_severe_chronic[i ,k]        <- df_nonhosp_severe[i, k]*lhs_sample$severe_chronic[k]
      df_nonhosp_severe_recov[i, k]          <- df_nonhosp_severe[i, k]*(1 - (lhs_sample$severe_chronic[k]))
      df_nonhosp_severe_chronic_mild[i, k]   <- df_nonhosp_severe[i, k]*lhs_sample$severe_mild[k]
      df_nonhosp_severe_chronic_severe[i, k] <- df_nonhosp_severe[i, k]*lhs_sample$severe_severe[k]
      
      df_nonhosp_mod[i, k]                   <- df_symp[i, k]*lhs_sample$nonhosp_mod[k]
      df_nonhosp_mod_chronic[i, k]           <- df_nonhosp_mod[i, k]*lhs_sample$mod_chronic[k]
      df_nonhosp_mod_recov[i, k]             <- df_nonhosp_mod[i, k]*(1 - (lhs_sample$mod_chronic[k]))
      df_nonhosp_mod_chronic_mild[i, k]      <- df_nonhosp_mod[i, k]*lhs_sample$mod_mild[k]
      df_nonhosp_mod_chronic_severe[i, k]    <- df_nonhosp_mod[i, k]*lhs_sample$mod_severe[k]
      
      df_nonhosp_mild[i, k]                   <- df_symp[i, k]*lhs_sample$nonhosp_mild[k]
      df_nonhosp_mild_chronic[i, k]           <- df_nonhosp_mild[i, k]*lhs_sample$mild_chronic[k]
      df_nonhosp_mild_recov[i, k]             <- df_nonhosp_mild[i, k]*(1 - (lhs_sample$mild_chronic[k]))
      df_nonhosp_mild_chronic_mild[i, k]      <- df_nonhosp_mild[i, k]
      
      df_hosp_fatal[i, k]                     <- df_hosp[i, k]*lhs_sample$fatal_hosp[k]
      
      # daly estimation
      
      df_yld_hosp_chronic[i, k]        <- (df_hosp_chronic_mild[i, k]*lhs_sample$dw_chronic_mild[k]*lhs_sample$dur_chronic[k]) + 
        (df_hosp_chronic_severe[i, k]*lhs_sample$dw_chronic_severe[k]*lhs_sample$dur_chronic[k])
      
      df_yld_hosp_acute[i, k]          <- df_hosp_recov[i, k]*lhs_sample$dw_hosp[k]*lhs_sample$dur_acute[k]
      
      df_yld_nh_severe_chronic[i, k]   <- (df_nonhosp_severe_chronic_mild[i, k]*lhs_sample$dw_chronic_mild[k]*lhs_sample$dur_chronic[k]) + 
        (df_nonhosp_severe_chronic_severe[i, k]*lhs_sample$dw_chronic_severe[k]*lhs_sample$dur_chronic[k])
      
      df_yld_nh_severe_acute[i, k]     <- df_nonhosp_severe_recov[i, k]*lhs_sample$dw_hosp[k]*lhs_sample$dur_acute[k]
      
      df_yld_nh_mod_chronic[i, k]      <- (df_nonhosp_mod_chronic_mild[i, k]*lhs_sample$dw_chronic_mild[k]*lhs_sample$dur_chronic[k]) + 
        (df_nonhosp_mod_chronic_severe[i, k]*lhs_sample$dw_chronic_severe[k]*lhs_sample$dur_chronic[k])
      
      df_yld_nh_mod_acute[i, k]        <- df_nonhosp_mod_recov[i, k]*lhs_sample$dw_hosp[k]*lhs_sample$dur_acute[k]
      
      df_yld_nh_mild_chronic[i, k]     <- (df_nonhosp_mild_chronic_mild[i, k]*lhs_sample$dw_chronic_mild[k]*lhs_sample$dur_chronic[k])
      
      df_yld_nh_mild_acute[i, k]        <- df_nonhosp_mild_recov[i, k]*lhs_sample$dw_hosp[k]*lhs_sample$dur_acute[k]
      
      # final daly 
      
      df_yld_acute[i, k]                <- df_yld_hosp_acute[i, k] + df_yld_nh_severe_acute[i, k] + df_yld_nh_mod_acute[i, k] + df_yld_nh_mild_acute[i, k] 
      df_yld_chronic[i, k]              <- df_yld_hosp_chronic[i, k] + df_yld_nh_severe_chronic[i, k] + df_yld_nh_mod_chronic[i, k] + df_yld_nh_mild_chronic[i, k]
      df_yll[i, k]                      <- df_hosp_fatal[i, k] * lhs_sample$le_lost[k]
      df_daly[i, k]                     <- df_yld_acute[i, k]  + df_yld_chronic[i, k]
      
    }
  }
  
  result_list <- list(symp = df_symp, asymp = df_asymp, df_hosp = df_hosp, df_hosp_recov = df_hosp_recov,
                      df_hosp_chronic = df_hosp_chronic, df_hosp_chronic_mild = df_hosp_chronic_mild,
                      df_hosp_chronic_severe = df_hosp_chronic_severe, df_nonhosp_severe = df_nonhosp_severe, df_nonhosp_severe_recov = df_nonhosp_severe_recov,
                      df_nonhosp_severe_chronic = df_nonhosp_severe_chronic, 
                      df_nonhosp_severe_chronic_mild = df_nonhosp_severe_chronic_mild, df_nonhosp_severe_chronic_severe = df_nonhosp_severe_chronic_severe,
                      df_nonhosp_mod = df_nonhosp_mod, df_nonhosp_mod_recov = df_nonhosp_mod_recov, df_nonhosp_mod_chronic = df_nonhosp_mod_chronic,
                      df_nonhosp_mod_chronic_mild = df_nonhosp_mod_chronic_mild, df_nonhosp_mod_chronic_severe = df_nonhosp_mod_chronic_severe,
                      df_nonhosp_mild = df_nonhosp_mild, df_nonhosp_mild_recov = df_nonhosp_mild_recov, df_nonhosp_mild_chronic = df_nonhosp_mild_chronic,
                      df_nonhosp_mild_chronic_mild = df_nonhosp_mild_chronic_mild, df_hosp_fatal = df_hosp_fatal,
                      df_yld_hosp_chronic = df_yld_hosp_chronic, df_yld_hosp_acute = df_yld_hosp_acute, 
                      df_yld_nh_severe_chronic = df_yld_nh_severe_chronic, df_yld_nh_severe_acute = df_yld_nh_severe_acute, 
                      df_yld_nh_mod_chronic = df_yld_nh_mod_chronic, df_yld_nh_mod_acute = df_yld_nh_mod_acute, 
                      df_yld_nh_mild_chronic = df_yld_nh_mild_chronic, df_yld_nh_mild_acute = df_yld_nh_mild_acute, 
                      df_yld_acute = df_yld_acute, df_yld_chronic = df_yld_chronic, df_yll = df_yll, df_daly = df_daly)  
  
  for (j in names(result_list)) {
    result_list[[j]] <- cbind(infection_df[, c( "country", "continent")], result_list[[j]])
    result_list[[j]]$mid <- apply(result_list[[j]][, c(3:1002)], 1, median, na.rm = TRUE)
    result_list[[j]]$lo  <- apply(result_list[[j]][, c(3:1002)], 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
    result_list[[j]]$hi  <- apply(result_list[[j]][, c(3:1002)], 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))
  }
  return(result_list)
}


## age group specific hospitalisation and fatality rate
age_specific_subinf_under40 <- function(infection_df, lhs_sample) {
  
  result_list <- list()  # Initialize a list to store results
  
  for (group in 1:4) {
    
    df_symp <- df_asymp <- df_hosp <- df_hosp_acute <- df_hosp_chronic_6m <- df_hosp_chronic_12m <- df_hosp_chronic_30m <- df_hosp_subac <- df_hosp_fatal <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
    colnames(df_symp) <- colnames(df_hosp) <- colnames(df_asymp) <- colnames(df_hosp_acute) <- colnames(df_hosp_chronic_6m) <- colnames(df_hosp_chronic_12m) <- colnames(df_hosp_chronic_30m) <- colnames(df_hosp_subac) <- colnames(df_hosp_fatal) <- paste0("Run_", 1:nrow(lhs_sample))
    
    df_nonhosp <- df_nonhosp_acute <- df_nonhosp_chronic_6m <- df_nonhosp_chronic_12m <- df_nonhosp_chronic_30m <- df_nonhosp_subac <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
    colnames(df_nonhosp) <- colnames(df_nonhosp_acute) <- colnames(df_nonhosp_chronic_6m) <- colnames(df_nonhosp_chronic_12m) <- colnames(df_nonhosp_chronic_30m) <- colnames(df_nonhosp_subac) <- paste0("Run_", 1:nrow(lhs_sample))
    
    df_yld_hosp_chronic <- df_yld_hosp_acute <- df_yld_hosp_subac <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
    colnames(df_yld_hosp_chronic) <- colnames(df_yld_hosp_acute) <- colnames(df_yld_hosp_subac) <- paste0("Run_", 1:nrow(lhs_sample))
    
    df_yld_nh_chronic <- df_yld_nh_acute <- df_yld_nh_subac <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
    colnames(df_yld_nh_chronic) <- colnames(df_yld_nh_acute) <- colnames(df_yld_nh_subac) <- paste0("Run_", 1:nrow(lhs_sample))
    
    df_yld_acute <- df_yld_subac <- df_yld_chronic <- df_yll <- df_daly <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
    colnames(df_yld_acute) <- colnames(df_yld_subac) <- colnames(df_yld_chronic) <- colnames(df_yll) <- colnames(df_daly) <- paste0("Run_", 1:nrow(lhs_sample))
    
    
    for(k in 1:nrow(lhs_sample)) {
      
      for(i in 1:nrow(infection_df)) {
        
        # Mid estimates
        med_inf <- infection_df[[paste0("tot_infec_med_group", group)]][i]
        
        # Apply Asian symptomatic to Asian continent 
        
        if (infection_df$continent[i] == "Asia") {
          
          df_symp[i, k] <- med_inf * lhs_sample$symp_asia[k]
          
          # Apply ECSA symptomatic to African continent
          
        } else if (infection_df$continent[i] == "Africa") {
          
          df_symp[i, k]  <- med_inf * lhs_sample$symp_africa[k]
          
          # Apply ECSA symptomatic to African continent
          
        } else if (infection_df$continent[i] == "North America") {
          
          df_symp[i, k]  <- med_inf * lhs_sample$symp_america[k]
          
        } else if (infection_df$continent[i] == "South America") {
          
          df_symp[i, k]  <- med_inf * lhs_sample$symp_america[k]
          
          # Apply overall symptomatic to other continents (EU/Oceania)
          
        } else {
          
          df_symp[i, k] <- med_inf * lhs_sample$symp_overall[k]
          
        }
        
        
        df_asymp[i, k]               <- med_inf - df_symp[i, k]
        #hosp_rate                    <- hosp_sample[[paste0("hosp_", group)]][k]
        df_hosp[i, k]                <- df_symp[i, k]*lhs_sample$hosp[k]
        df_hosp_acute[i, k]          <- df_hosp[i, k]*lhs_sample$acute[k]
        df_hosp_chronic_6m[i, k]     <- df_hosp[i, k]*lhs_sample$chr6m[k] 
        df_hosp_chronic_12m[i, k]    <- df_hosp[i, k]*lhs_sample$chr12m[k] 
        df_hosp_chronic_30m[i, k]    <- df_hosp[i, k]*lhs_sample$chr30m[k]
        df_hosp_subac[i, k]          <- df_hosp[i, k]*lhs_sample$subac[k]
        
        df_nonhosp[i, k]              <- df_symp[i, k] - df_hosp[i, k]
        df_nonhosp_acute[i, k]        <- df_nonhosp[i, k]*lhs_sample$acute[k]
        df_nonhosp_chronic_6m[i ,k]   <- df_nonhosp[i, k]*lhs_sample$chr6m[k]
        df_nonhosp_chronic_12m[i ,k]  <- df_nonhosp[i, k]*lhs_sample$chr12m[k] 
        df_nonhosp_chronic_30m[i ,k]  <- df_nonhosp[i, k]*lhs_sample$chr30m[k]
        df_nonhosp_subac[i, k]        <- df_nonhosp[i, k]*lhs_sample$subac[k]
        
        fatal_rate_hosp               <- fatal_sample[[paste0("fatal_", group)]][k]
        fatal_rate_nonhosp            <- nh_fatal_sample[[paste0("nh_fatal_", group)]][k]
        df_hosp_fatal[i, k]           <- df_hosp[i, k]*fatal_rate_hosp + df_nonhosp[i, k]*fatal_rate_nonhosp
        #df_hosp_fatal[i, k]          <- df_hosp[i, k]*lhs_sample$fatal_hosp[k] + df_nonhosp[i, k]*lhs_sample$fatal_nonhosp[k]
        
        # YLD estimation
        
        df_yld_hosp_chronic[i, k]        <- df_hosp_chronic_6m[i, k]*lhs_sample$dur_6m[k]*lhs_sample$dw_chronic[k] + df_hosp_chronic_12m[i, k]*lhs_sample$dur_12m[k]*lhs_sample$dw_chronic[k] + df_hosp_chronic_30m[i, k]*lhs_sample$dur_30m[k]*lhs_sample$dw_chronic[k]
        df_yld_hosp_acute[i, k]          <- df_hosp_acute[i, k]*lhs_sample$dw_hosp[k]*lhs_sample$dur_acute[k]
        df_yld_hosp_subac[i, k]          <- df_hosp_subac[i, k]*lhs_sample$dw_subac[k]*lhs_sample$dur_subac[k]
        
        df_yld_nh_chronic[i, k]      <- df_nonhosp_chronic_6m[i, k]*lhs_sample$dur_6m[k]*lhs_sample$dw_chronic[k] + df_nonhosp_chronic_12m[i, k]*lhs_sample$dur_12m[k]*lhs_sample$dw_chronic[k] + df_nonhosp_chronic_30m[i, k]*lhs_sample$dur_30m[k]*lhs_sample$dw_chronic[k]
        df_yld_nh_acute[i, k]        <- df_nonhosp_acute[i, k]*lhs_sample$dw_hosp[k]*lhs_sample$dur_acute[k]
        df_yld_nh_subac[i, k]        <- df_nonhosp_subac[i, k]*lhs_sample$dw_subac[k]*lhs_sample$dur_subac[k]
        
        # final daly 
        
        df_yld_acute[i, k]                <- df_yld_hosp_acute[i, k] + df_yld_nh_acute[i, k] 
        df_yld_subac[i, k]                <- df_yld_hosp_subac[i, k] + df_yld_nh_subac[i, k] 
        df_yld_chronic[i, k]              <- df_yld_hosp_chronic[i, k] + df_yld_nh_chronic[i, k] 
        
        le_left                           <- le_sample[[paste0("le_", group)]][k]
        df_yll[i, k]                      <- df_hosp_fatal[i, k]*le_left
        df_daly[i, k]                     <- df_yld_acute[i, k] + df_yld_subac[i, k] + df_yld_chronic[i, k] + df_yll[i, k] 
        
      }
    }
    
    result_list[[paste0("Group_", group)]] <- list(
      symp = df_symp, 
      asymp = df_asymp, 
      df_hosp = df_hosp, 
      df_hosp_acute = df_hosp_acute,
      df_hosp_chronic_6m = df_hosp_chronic_6m,
      df_hosp_chronic_12m = df_hosp_chronic_12m,
      df_hosp_chronic_30m = df_hosp_chronic_30m, 
      df_hosp_subac = df_hosp_subac, 
      df_hosp_fatal = df_hosp_fatal,
      df_nonhosp = df_nonhosp, 
      df_nonhosp_acute = df_nonhosp_acute,
      df_nonhosp_chronic_6m = df_nonhosp_chronic_6m, 
      df_nonhosp_chronic_12m = df_nonhosp_chronic_12m,
      df_nonhosp_chronic_30m = df_nonhosp_chronic_30m,
      df_nonhosp_subac = df_nonhosp_subac,
      df_yld_hosp_chronic = df_yld_hosp_chronic, 
      df_yld_hosp_acute = df_yld_hosp_acute, 
      df_yld_hosp_subac = df_yld_hosp_subac,
      df_yld_nh_chronic = df_yld_nh_chronic, 
      df_yld_nh_acute = df_yld_nh_acute, 
      df_yld_nh_subac = df_yld_nh_subac,
      df_yld_acute = df_yld_acute, 
      df_yld_subac = df_yld_subac, 
      df_yld_chronic = df_yld_chronic, 
      df_yll = df_yll, 
      df_daly = df_daly
    )
    
    # Add country and continent columns, and compute statistics for each dataset
    for (j in names(result_list[[paste0("Group_", group)]])) {
      result_list[[paste0("Group_", group)]][[j]] <- cbind(infection_df[, c("country", "continent", "tot_pop")], result_list[[paste0("Group_", group)]][[j]])
      result_list[[paste0("Group_", group)]][[j]]$mid <- apply(result_list[[paste0("Group_", group)]][[j]][, 4:(3 + nrow(lhs_sample))], 1, median, na.rm = TRUE)
      result_list[[paste0("Group_", group)]][[j]]$lo <- apply(result_list[[paste0("Group_", group)]][[j]][, 4:(3 + nrow(lhs_sample))], 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
      result_list[[paste0("Group_", group)]][[j]]$hi <- apply(result_list[[paste0("Group_", group)]][[j]][, 4:(3 + nrow(lhs_sample))], 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))
    }
    
  }
  
  return(result_list)
}
age_specific_subinf_over40 <- function(infection_df, lhs_sample) {
  
  result_list <- list()  # Initialize a list to store results
  
  for (group in 5:9) {
    
    df_symp <- df_asymp <- df_hosp <- df_hosp_acute <- df_hosp_chronic_6m <- df_hosp_chronic_12m <- df_hosp_chronic_30m <- df_hosp_subac <- df_hosp_fatal <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
    colnames(df_symp) <- colnames(df_hosp) <- colnames(df_asymp) <- colnames(df_hosp_acute) <- colnames(df_hosp_chronic_6m) <- colnames(df_hosp_chronic_12m) <- colnames(df_hosp_chronic_30m) <- colnames(df_hosp_subac) <- colnames(df_hosp_fatal) <- paste0("Run_", 1:nrow(lhs_sample))
    
    df_nonhosp <- df_nonhosp_acute <- df_nonhosp_chronic_6m <- df_nonhosp_chronic_12m <- df_nonhosp_chronic_30m <- df_nonhosp_subac <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
    colnames(df_nonhosp) <- colnames(df_nonhosp_acute) <- colnames(df_nonhosp_chronic_6m) <- colnames(df_nonhosp_chronic_12m) <- colnames(df_nonhosp_chronic_30m) <- colnames(df_nonhosp_subac) <- paste0("Run_", 1:nrow(lhs_sample))
    
    df_yld_hosp_chronic <- df_yld_hosp_acute <- df_yld_hosp_subac <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
    colnames(df_yld_hosp_chronic) <- colnames(df_yld_hosp_acute) <- colnames(df_yld_hosp_subac) <- paste0("Run_", 1:nrow(lhs_sample))
    
    df_yld_nh_chronic <- df_yld_nh_acute <- df_yld_nh_subac <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
    colnames(df_yld_nh_chronic) <- colnames(df_yld_nh_acute) <- colnames(df_yld_nh_subac) <- paste0("Run_", 1:nrow(lhs_sample))
    
    df_yld_acute <- df_yld_subac <- df_yld_chronic <- df_yll <- df_daly <- data.frame(matrix(NA, nrow = nrow(infection_df), ncol = nrow(lhs_sample)))
    colnames(df_yld_acute) <- colnames(df_yld_subac) <- colnames(df_yld_chronic) <- colnames(df_yll) <- colnames(df_daly) <- paste0("Run_", 1:nrow(lhs_sample))
    
    
    for(k in 1:nrow(lhs_sample)) {
      
      for(i in 1:nrow(infection_df)) {
        
        # Mid estimates
        med_inf <- infection_df[[paste0("tot_infec_med_group", group)]][i]
        
        # Apply Asian symptomatic to Asian continent 
        
        if (infection_df$continent[i] == "Asia") {
          
          df_symp[i, k] <- med_inf * lhs_sample$symp_asia[k]
          
          # Apply ECSA symptomatic to African continent
          
        } else if (infection_df$continent[i] == "Africa") {
          
          df_symp[i, k]  <- med_inf * lhs_sample$symp_africa[k]
          
          # Apply ECSA symptomatic to African continent
          
        } else if (infection_df$continent[i] == "North America") {
          
          df_symp[i, k]  <- med_inf * lhs_sample$symp_america[k]
          
        } else if (infection_df$continent[i] == "South America") {
          
          df_symp[i, k]  <- med_inf * lhs_sample$symp_america[k]
          
          # Apply overall symptomatic to other continents (EU/Oceania)
          
        } else {
          
          df_symp[i, k] <- med_inf * lhs_sample$symp_overall[k]
          
        }
        
        
        df_asymp[i, k]               <- med_inf - df_symp[i, k]
        #hosp_rate                    <- hosp_sample[[paste0("hosp_", group)]][k]
        df_hosp[i, k]                <- df_symp[i, k]*lhs_sample$hosp[k]
        df_hosp_acute[i, k]          <- df_hosp[i, k]*lhs_sample$acute[k]
        df_hosp_chronic_6m[i, k]     <- df_hosp[i, k]*lhs_sample$chr6m[k] 
        df_hosp_chronic_12m[i, k]    <- df_hosp[i, k]*lhs_sample$chr12m[k] 
        df_hosp_chronic_30m[i, k]    <- df_hosp[i, k]*lhs_sample$chr30m[k]
        df_hosp_subac[i, k]          <- df_hosp[i, k]*lhs_sample$subac[k]
        
        #fatal_rate                   <- fatal_sample[[paste0("fatal_", group)]][k]

        df_nonhosp[i, k]              <- df_symp[i, k] - df_hosp[i, k]
        df_nonhosp_acute[i, k]        <- df_nonhosp[i, k]*lhs_sample$acute[k]
        df_nonhosp_chronic_6m[i ,k]   <- df_nonhosp[i, k]*lhs_sample$chr6m[k]
        df_nonhosp_chronic_12m[i ,k]  <- df_nonhosp[i, k]*lhs_sample$chr12m[k] 
        df_nonhosp_chronic_30m[i ,k]  <- df_nonhosp[i, k]*lhs_sample$chr30m[k]
        df_nonhosp_subac[i, k]        <- df_nonhosp[i, k]*lhs_sample$subac[k]
        
        fatal_rate_hosp               <- fatal_sample[[paste0("fatal_", group)]][k]
        fatal_rate_nonhosp            <- nh_fatal_sample[[paste0("nh_fatal_", group)]][k]
        df_hosp_fatal[i, k]           <- df_hosp[i, k]*fatal_rate_hosp + df_nonhosp[i, k]*fatal_rate_nonhosp
        
        #df_hosp_fatal[i, k]          <- df_hosp[i, k]*lhs_sample$fatal_hosp[k] + df_nonhosp[i, k]*lhs_sample$fatal_nonhosp[k]
        
        # YLD estimation
        
        df_yld_hosp_chronic[i, k]        <- df_hosp_chronic_6m[i, k]*lhs_sample$dur_6m[k]*lhs_sample$dw_chronic[k] + df_hosp_chronic_12m[i, k]*lhs_sample$dur_12m[k]*lhs_sample$dw_chronic[k] + df_hosp_chronic_30m[i, k]*lhs_sample$dur_30m[k]*lhs_sample$dw_chronic[k]
        df_yld_hosp_acute[i, k]          <- df_hosp_acute[i, k]*lhs_sample$dw_hosp[k]*lhs_sample$dur_acute[k]
        df_yld_hosp_subac[i, k]          <- df_hosp_subac[i, k]*lhs_sample$dw_subac[k]*lhs_sample$dur_subac[k]
        
        df_yld_nh_chronic[i, k]      <- df_nonhosp_chronic_6m[i, k]*lhs_sample$dur_6m[k]*lhs_sample$dw_chronic[k] + df_nonhosp_chronic_12m[i, k]*lhs_sample$dur_12m[k]*lhs_sample$dw_chronic[k] + df_nonhosp_chronic_30m[i, k]*lhs_sample$dur_30m[k]*lhs_sample$dw_chronic[k]
        df_yld_nh_acute[i, k]        <- df_nonhosp_acute[i, k]*lhs_sample$dw_hosp[k]*lhs_sample$dur_acute[k]
        df_yld_nh_subac[i, k]        <- df_nonhosp_subac[i, k]*lhs_sample$dw_subac[k]*lhs_sample$dur_subac[k]
        
        # final daly 
        
        df_yld_acute[i, k]                <- df_yld_hosp_acute[i, k] + df_yld_nh_acute[i, k] 
        df_yld_subac[i, k]                <- df_yld_hosp_subac[i, k] + df_yld_nh_subac[i, k] 
        df_yld_chronic[i, k]              <- df_yld_hosp_chronic[i, k] + df_yld_nh_chronic[i, k] 
        
        le_left                           <- le_sample[[paste0("le_", group)]][k]
        df_yll[i, k]                      <- df_hosp_fatal[i, k]*le_left  
        
        df_daly[i, k]                     <- df_yld_acute[i, k] + df_yld_subac[i, k] + df_yld_chronic[i, k] + df_yll[i, k] 
        
      }
    }
    
    result_list[[paste0("Group_", group)]] <- list(
      symp = df_symp, 
      asymp = df_asymp, 
      df_hosp = df_hosp, 
      df_hosp_acute = df_hosp_acute,
      df_hosp_chronic_6m = df_hosp_chronic_6m,
      df_hosp_chronic_12m = df_hosp_chronic_12m,
      df_hosp_chronic_30m = df_hosp_chronic_30m, 
      df_hosp_subac = df_hosp_subac, 
      df_hosp_fatal = df_hosp_fatal,
      df_nonhosp = df_nonhosp, 
      df_nonhosp_acute = df_nonhosp_acute,
      df_nonhosp_chronic_6m = df_nonhosp_chronic_6m, 
      df_nonhosp_chronic_12m = df_nonhosp_chronic_12m,
      df_nonhosp_chronic_30m = df_nonhosp_chronic_30m,
      df_nonhosp_subac = df_nonhosp_subac,
      df_yld_hosp_chronic = df_yld_hosp_chronic, 
      df_yld_hosp_acute = df_yld_hosp_acute, 
      df_yld_hosp_subac = df_yld_hosp_subac,
      df_yld_nh_chronic = df_yld_nh_chronic, 
      df_yld_nh_acute = df_yld_nh_acute, 
      df_yld_nh_subac = df_yld_nh_subac,
      df_yld_acute = df_yld_acute, 
      df_yld_subac = df_yld_subac, 
      df_yld_chronic = df_yld_chronic, 
      df_yll = df_yll, 
      df_daly = df_daly
    )
    
    # Add country and continent columns, and compute statistics for each dataset
    for (j in names(result_list[[paste0("Group_", group)]])) {
      result_list[[paste0("Group_", group)]][[j]] <- cbind(infection_df[, c("country", "continent", "tot_pop")], result_list[[paste0("Group_", group)]][[j]])
      result_list[[paste0("Group_", group)]][[j]]$mid <- apply(result_list[[paste0("Group_", group)]][[j]][, 4:(3 + nrow(lhs_sample))], 1, median, na.rm = TRUE)
      result_list[[paste0("Group_", group)]][[j]]$lo <- apply(result_list[[paste0("Group_", group)]][[j]][, 4:(3 + nrow(lhs_sample))], 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
      result_list[[paste0("Group_", group)]][[j]]$hi <- apply(result_list[[paste0("Group_", group)]][[j]][, 4:(3 + nrow(lhs_sample))], 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))
    }
    
  }
  
  return(result_list)
}


# post processing sub burdens
sub_inf_global_count <- function(psa_df){
  
  global <- psa_df %>%
    ungroup() %>%
    summarise(
      tot_med  = sum(mid, na.rm = TRUE),
      tot_lo   = sum(lo, na.rm = TRUE),
      tot_hi   = sum(hi, na.rm = TRUE),
      tot_pop  = sum(tot_pop, na.rm = TRUE)
    ) %>% 
    as.data.frame()
  
  cont <- psa_df %>% group_by(continent) %>%
    summarise(
      tot_med  = sum(mid, na.rm = TRUE),
      tot_lo   = sum(lo, na.rm = TRUE),
      tot_hi   = sum(hi, na.rm = TRUE),
      continent      = first(continent),
      tot_pop  = sum(tot_pop, na.rm = TRUE)
    ) %>% 
    as.data.frame()
  
  country <- psa_df %>% group_by(country) %>%
    summarise(
      tot_med  = sum(mid, na.rm = TRUE),
      tot_lo   = sum(lo, na.rm = TRUE),
      tot_hi   = sum(hi, na.rm = TRUE),
      continent      = first(continent),
      tot_pop  = sum(tot_pop, na.rm = TRUE)
    ) %>% 
    as.data.frame()
  
  return(list(global  = global,
              cont    = cont,
              country = country))
  
}

postprocess_sub_burden <- function(sub_burden_psa, infection_df, sub_inf_global_count) {
  
  # Add total pop
  
  for(category in names(sub_burden_psa)) {
    
    sub_burden_psa[[category]]$tot_pop <- infection_df$tot_pop
  }
  
  # Apply sub_inf_glob_count function
  global_count <- lapply(sub_burden_psa, sub_inf_global_count)
  
  # Extract global, continental, country results
  
  global_results <- list()
  cont_results <- list()
  country_results <- list()
  
  for(category in names(global_count)) {
    
    global_results[[category]] <- global_count[[category]]$global
    cont_results[[category]] <- global_count[[category]]$cont
    country_results[[category]] <- global_count[[category]]$country
    
  }
  
  # Compute rates per million for continental data
  #for (category in names(cont_results)) {
  #  cont_results[[category]]$rate_mid <- (cont_results[[category]]$tot_med / cont_results[[category]]$tot_pop) * 1000000
  #  cont_results[[category]]$rate_lo <- (cont_results[[category]]$tot_lo / cont_results[[category]]$tot_pop) * 1000000
  #  cont_results[[category]]$rate_hi <- (cont_results[[category]]$tot_hi / cont_results[[category]]$tot_pop) * 1000000
  #}
  
  # Add type labels to global and continental 
  
  type <- c(
    symp                   = "symptomatic", 
    asymp                  = "asymptomatic", 
    df_hosp                = "hospitalisation", 
    df_hosp_acute          = "hospitalisation_acute",
    df_hosp_chronic_6m     = "hospitalisation_chronic6m",
    df_hosp_chronic_12m    = "hospitalisation_chronic12m",
    df_hosp_chronic_30m    = "hospitalisation_chronic30m", 
    df_hosp_subac          = "hospitalisation_subacute", 
    df_hosp_fatal          = "fatal",
    df_nonhosp             = "nonhospitalisation", 
    df_nonhosp_acute       = "nonhosp_acute",
    df_nonhosp_chronic_6m  = "nonhosp_chronic_6m", 
    df_nonhosp_chronic_12m = "nonhosp_chronic_12m",
    df_nonhosp_chronic_30m = "nonhosp_chronic_30m",
    df_nonhosp_subac       = "nonhosp_subac",
    df_yld_hosp_chronic    = "yld_hosp_chronic", 
    df_yld_hosp_acute      = "yld_hosp_acute", 
    df_yld_hosp_subac      = "yld_hosp_subac",
    df_yld_nh_chronic      = "yld_nh_chronic", 
    df_yld_nh_acute        = "yld_nh_acute", 
    df_yld_nh_subac        = "yld_nh_subac",
    df_yld_acute           = "yld_acute", 
    df_yld_subac           = "yld_subac", 
    df_yld_chronic         = "yld_chronic", 
    df_yll                 = "yll", 
    df_daly                = "daly"
  )
  
  
  for(category in names(global_results)) {
    
    global_results[[category]]$type <- type[[category]]
    
  }
  
  for(category in names(cont_results)) {
    
    cont_results[[category]]$type <- type[[category]]
    
  }
  
  for(category in names(cont_results)) {
    
    cont_results[[category]]$type <- type[[category]]
    
  }
  
  for(category in names(country_results)) {
    
    country_results[[category]]$type <- type[[category]]
    
  }
  
  
  # Append all global and continental results 
  
  all_results_global    <- do.call(rbind, global_results)
  all_results_continent <- do.call(rbind, cont_results)
  all_results_country   <- do.call(rbind, country_results)
  
  return(list(global_count = global_count, global = global_results, 
              cont = cont_results, country = country_results, 
              all_results_global = all_results_global, all_results_continent = all_results_continent,
              all_results_country = all_results_country))
  
}

