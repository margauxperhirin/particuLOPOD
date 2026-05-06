#' =============================================================================
#' @name model_continuous
#' @description sub-pipeline corresponding to the model fitting procedure for
#' continuous data. Called via the model_wrapper.R function, thus no default
#' parameter values are defined.
#' 
#' @param CALL the call object from the master pipeline. 
#' @param QUERY the query object from the master pipeline.
#' 
#' @return returns a list object containing the best model, associated hyper
#' parameters and predicted values per re sampling folds

model_continuous <- function(CALL, QUERY){
  
  # --- 1. Initialize
  # --- 1.1. Store hyperparameter set in MODEL object
  MODEL <- CALL$HP
  
  # --- 1.2. Get model list
  model_name <- CALL$MODEL_LIST
  rm(CALL) # Remove CALL for memory use
  gc()
  
  targets <- c("psd1","psd2","psd3")
  
  # --- 2. Loop with all selected models
  for(i in seq_along(model_name)){
    
    message(paste(Sys.time(), "--- Start hyper parameter tuning for", model_name[i], "---"))
    
    if(model_name[i] == "RF"){
      
      library(randomForestSRC)
      
      grid <- MODEL[[model_name[i]]][["model_grid"]]
      resamples <- QUERY$FOLDS$resample_split
      targets <- c("psd1", "psd2", "psd3")
      
      # --- Data full
      full_df <- cbind(QUERY$Y[, targets], QUERY$X[, QUERY$SUBFOLDER_INFO$ENV_VAR])
      
      # 1. Préparation des Folds (Calculé une seule fois)
      message("--- Pre-processing folds for randomForestSRC ---")
      prepped_folds <- lapply(seq_len(nrow(resamples)), function(fold) {
        split <- resamples$splits[[fold]]
        train_data <- analysis(split)
        test_data  <- assessment(split)
        
        list(train_df = train_data[, c(targets, QUERY$SUBFOLDER_INFO$ENV_VAR)],
             test_df  = test_data[, QUERY$SUBFOLDER_INFO$ENV_VAR],
             truth    = as.matrix(test_data[, targets]))
      })
      
      # 2. Création dynamique de la formule Multivariée
      rf_formula <- as.formula(paste0("Multivar(", paste(targets, collapse = ", "), ") ~ ."))
      
      results <- list()
      
      message(paste(Sys.time(), "--- Starting Sequential Grid Search (randomForestSRC) ---"))
      
      # 3. Grid search classique (Boucle For)
      for (g in seq_len(nrow(grid))) {
        
        params <- grid[g, ]
        fold_rmse <- numeric(length(prepped_folds))
        
        for (fold in seq_along(prepped_folds)) {
          fd <- prepped_folds[[fold]]
          
          # Fit du modèle Multivarié
          mod <- randomForestSRC::rfsrc(formula = rf_formula,
                                        data = fd$train_df,
                                        ntree = params$trees,
                                        nodesize = params$min_n,
                                        importance = "none") # Désactivé ici pour accélérer le grid search
          
          # Prédiction
          pred_obj <- predict(mod, newdata = fd$test_df)
          
          # Reconstruction de la matrice de prédiction
          pred_matrix <- sapply(targets, function(t) {
            pred_obj$regrOutput[[t]]$predicted
          })
          
          # Calcul RMSE
          fold_rmse[fold] <- sqrt(mean((fd$truth - pred_matrix)^2))
        }
        
        results[[g]] <- list(
          params = params,
          rmse = mean(fold_rmse)
        )
        
        message(paste("Grid", g, "/", nrow(grid), "- RMSE:", round(results[[g]]$rmse, 4)))
      }
      
      # 4. Sélection des meilleurs paramètres
      rmse_vals <- sapply(results, function(x) x$rmse)
      best_id <- which.min(rmse_vals)
      best_model <- results[[best_id]]
      
      message(paste("Best Grid - RMSE:", round(best_model$rmse, 4)))
      
      # 5. Fit Final (sur toutes les données)
      final_rf_obj <- randomForestSRC::rfsrc(formula = rf_formula,
                                             data = full_df,
                                             ntree = best_model$params$trees,
                                             nodesize = best_model$params$min_n,
                                             importance = TRUE) # Activé pour le modèle final
      
      
      # Prédictions finales
      final_pred_obj <- predict(final_rf_obj, newdata = full_df)
      final_pred_matrix <- sapply(targets, function(t) {
        final_pred_obj$regrOutput[[t]]$predicted
      })
      
      # 6. Stockage des résultats
      MODEL[[model_name[i]]][["best_params"]] <- best_model$params
      MODEL[[model_name[i]]][["best_rmse"]] <- best_model$rmse
      MODEL[[model_name[i]]][["final_fit"]] <- final_pred_matrix 
      MODEL[[model_name[i]]][["rfsrc_model"]] <- final_rf_obj 
      
    } else if (model_name[i] == "MLP") {
      
      library(nnet)
      
      grid <- MODEL[["MLP"]][["model_grid"]]
      resamples <- QUERY$FOLDS$resample_split
      targets <- c("psd1", "psd2", "psd3")
      
      # --- Data full 
      full_df <- cbind(QUERY$Y[, targets, drop = FALSE], QUERY$X_norm[, QUERY$SUBFOLDER_INFO$ENV_VAR, drop = FALSE])
      
      # 1. Préparation des Folds
      message("--- Pre-processing folds for nnet (MLP) ---")
      prepped_folds <- lapply(seq_len(nrow(resamples)), function(fold) {
        split <- resamples$splits[[fold]]
        train_data <- analysis(split)
        test_data  <- assessment(split)
        
        for(v in QUERY$SUBFOLDER_INFO$ENV_VAR) {
          train_data[[v]] <- (train_data[[v]] - QUERY$X_norm_params$means[v]) / QUERY$X_norm_params$sds[v]
          test_data[[v]]  <- (test_data[[v]]  - QUERY$X_norm_params$means[v]) / QUERY$X_norm_params$sds[v]
        }
        
        list(train_df = train_data[, c(targets, QUERY$SUBFOLDER_INFO$ENV_VAR)],
             test_df  = test_data[, QUERY$SUBFOLDER_INFO$ENV_VAR],
             truth    = as.matrix(test_data[, targets]))
      })
      
      nnet_formula <- as.formula(paste0("cbind(", paste(targets, collapse = ", "), ") ~ ."))
      
      results <- list()
      message(paste(Sys.time(), "--- Starting Sequential Grid Search (nnet) ---"))
      
      # 3. Grid search classique AVEC sauvegarde des prédictions
      for (g in seq_len(nrow(grid))) {
        
        params <- grid[g, ]
        fold_rmse <- numeric(length(prepped_folds))
        
        # Listes temporaires pour stocker les preds/obs de chaque fold
        fold_preds <- list()
        fold_obs <- list()
        
        for (fold in seq_along(prepped_folds)) {
          fd <- prepped_folds[[fold]]
          
          mod <- nnet::nnet(formula = nnet_formula,
                            data = fd$train_df,
                            size = params$hidden_units,
                            linout = TRUE,
                            trace = FALSE,
                            maxit = 1000)
          
          pred_matrix <- as.matrix(predict(mod, newdata = fd$test_df))
          
          fold_rmse[fold] <- sqrt(mean((fd$truth - pred_matrix)^2))
          fold_preds[[fold]] <- pred_matrix
          fold_obs[[fold]] <- fd$truth
        }
        
        # On sauvegarde les preds et obs combinées pour cette configuration de grille
        results[[g]] <- list(params = params, 
                             rmse = mean(fold_rmse),
                             preds = do.call(rbind, fold_preds),
                             obs = do.call(rbind, fold_obs))
        
        message(paste("Grid", g, "/", nrow(grid), "- RMSE:", round(results[[g]]$rmse, 4)))
      }
      
      # 4. Sélection des meilleurs paramètres
      rmse_vals <- sapply(results, function(x) x$rmse)
      best_id <- which.min(rmse_vals)
      best_model <- results[[best_id]]
      
      message(paste("Best Grid - RMSE:", round(best_model$rmse, 4)))
      
      # 5. Fit Final (sur toutes les données, pour utilisation future)
      final_nnet_obj <- nnet::nnet(formula = nnet_formula,
                                   data = full_df,
                                   size = best_model$params$hidden_units,
                                   linout = TRUE,
                                   trace = FALSE,
                                   maxit = 1000)
      
      final_pred_matrix <- predict(final_nnet_obj, newdata = full_df)
      colnames(final_pred_matrix) <- targets
      
      # 6. Stockage des résultats (Ajout des OOF preds/obs)
      MODEL[[model_name[i]]][["best_params"]] <- best_model$params
      MODEL[[model_name[i]]][["best_rmse"]]   <- best_model$rmse
      MODEL[[model_name[i]]][["final_fit"]]   <- final_pred_matrix 
      MODEL[[model_name[i]]][["nnet_model"]]  <- final_nnet_obj 
      MODEL[[model_name[i]]][["oof_preds"]]   <- best_model$preds
      MODEL[[model_name[i]]][["oof_obs"]]     <- best_model$obs
      
    } else {
      
      # --- 2.1. Define the formula
    # General formula for most models; special handling for GAM
    if(model_name[i] != "GAM"){
      formula_vars <- QUERY$SUBFOLDER_INFO$ENV_VAR %>% paste(collapse = " + ")
    } else {
      # GAM requires spline terms
      formula_vars <- paste0("s(", QUERY$SUBFOLDER_INFO$ENV_VAR %>% paste(collapse = ", k = 3) + s("), ", k = 3)")
    }
    formula <- paste0("measurementvalue ~ ", formula_vars) %>% as.formula()
    
    # --- 2.2. Define workflow adapted to hyper parameter tuning
    model_wf <- workflow() %>% 
      add_variables(outcomes = "measurementvalue", predictors = QUERY$SUBFOLDER_INFO$ENV_VAR) %>% 
      add_model(MODEL[[model_name[i]]][["model_spec"]], formula = formula)
    
    # --- 2.3. Run the model for each fold x (hyper parameter grid rows)
    # Runs hyper parameter tuning if a grid is present in the HP (e.g. no GLM tune)
    if(!is.null(MODEL[[model_name[i]]][["model_grid"]])){
      model_res <- model_wf %>% 
        tune_grid(resamples = QUERY$FOLDS$resample_split,
                  grid = MODEL[[model_name[i]]][["model_grid"]],
                  metrics = yardstick::metric_set(rmse),
                  control = control_grid(verbose = TRUE, allow_par = FALSE))
    } else {
      model_res <- model_wf %>% 
        fit_resamples(resamples = QUERY$FOLDS$resample_split,
                      control = control_resamples(verbose = TRUE, allow_par = FALSE))
    }
    
    # --- 2.4. Select best hyper parameter set
    # Based on RMSE values per model run (rsq does not work with 0's)
    model_best <- model_res %>% select_best(metric = "rmse")
    # Retrieve the corresponding RMSE as well 
    MODEL[[model_name[i]]][["best_fit"]] <- model_res %>% show_best(metric = "rmse") %>% .[1,]
    
    # --- 2.5. Define final workflow
    final_wf <- model_wf %>% finalize_workflow(model_best)
    
    # --- 2.6. Run the model on same cross validation splits
    # We have one fit per cross validation saved in a list, to be passed to further steps
    final_fit <- lapply(seq_along(QUERY$FOLDS$resample_folds), function(x){
      out <- final_wf %>% last_fit(QUERY$FOLDS$resample_split$splits[[x]])
      return(out)
    }) # loop over the same cross validation folds
    
    MODEL[[model_name[i]]][["final_wf"]] <- final_wf
    MODEL[[model_name[i]]][["final_fit"]] <- final_fit
    
    }
    
    # --- 2.7. Display information
    message(paste(Sys.time(), "--- DONE ---"))
  }
  
  return(MODEL)
  
} # END FUNCTION
