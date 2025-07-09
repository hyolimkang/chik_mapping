default_rmse <- list()

for (i in 1:length(merge_list)) {
  merge_data <- as.data.frame(merge_list[[i]])
  
  foi_rf <- randomForest(
    logfoi ~ Tsuit + PRCP + GDP + Albo + Aegyp + CHIKRisk, 
    data = merge_data,
    mtry = floor(n_features / 3),
    importance = TRUE,
    seed = 123
  )
  
  log_oob_predictions <- predict(foi_rf, merge_data)
  
  result <- sqrt(mean((merge_data$logfoi - log_oob_predictions)^2))
  
  default_rmse[[i]] <- result
}


model_results <- list()

for (j in 1:length(train_list)) {
  train_data <- as.data.frame(train_list[[j]])
  
  best_rmse <- default_rmse[[j]]
  best_params <- NULL 
  
  for (i in seq_len(nrow(hyper_grid))) {
    # Fit model for ith hyperparameter combination using randomForest
    fit <- randomForest(
      logfoi ~ Tsuit + PRCP + GDP + Albo + Aegyp + CHIKRisk,
      data = train_data,
      ntree = n_features * 10,  # equivalent to num.trees in ranger
      mtry = hyper_grid$mtry[i],
      nodesize = hyper_grid$min.node.size[i],  # equivalent to min.node.size in ranger
      replace = hyper_grid$replace[i],
      sampsize = ifelse(hyper_grid$replace[i], nrow(train_data), nrow(train_data) * hyper_grid$sample.fraction[i]),  # handling sample size
      keep.inbag = TRUE,
      importance = TRUE,
      do.trace = FALSE,
      seed = 123  # setting seed
    )
    
    # current prediction error of ith model
    log_oob_pred <- predict(fit, train_data)
    current_rmse <- sqrt(mean((train_data$logfoi - log_oob_pred)^2))
    
    hyper_grid$rmse[i] <- current_rmse 
    hyper_grid$perc_gain[i] <- (default_rmse[[j]] - current_rmse) / default_rmse[[j]] * 100
  }
  
  best_params_df <- hyper_grid %>%
    dplyr::arrange(desc(perc_gain)) %>%
    dplyr::slice_head(n = 1)
  
  # Extracting best parameters
  best_params <- best_params_df[1, ]
  best_rmse <- best_params_df$rmse
  
  # Store the best parameters and RMSE for each model
  model_results[[j]] <- list(best_params = best_params, best_rmse = best_rmse)
}