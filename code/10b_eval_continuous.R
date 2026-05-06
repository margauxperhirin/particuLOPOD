#' =============================================================================
#' @name eval_continuous
#' @description sub-pipeline for model evaluation corresponding to continuous data
#' @param CALL the call object from the master pipeline
#' @param QUERY the query object from the master pipeline
#' @param MODEL the models object from the master pipeline
#' @param PDF_PATH preset path for the pdf save
#' @return the MODEL object updated with evaluation metric values (model performance
#' metric and variable importance metric)
#' @return variable importance plots as PDF file

eval_continuous <- function(CALL,
                            QUERY,
                            MODEL,
                            PDF_PATH){
    
    # --- 1. Model performance assessment
    for(i in MODEL$MODEL_LIST){
      
      if(i == "RF"){
        # --- 1.1 SPECIFIQUE RF : R2 basée sur les prédictions OOB 
        rf_mod <- MODEL[[i]][["rfsrc_model"]]
        
        # Extraction des prédictions Out-Of-Bag pour chaque cible
        pred_list <- lapply(c("psd1", "psd2", "psd3"), function(t) {
          rf_mod$regrOutput[[t]]$predicted.oob
        })
        
        pred <- do.call(cbind, pred_list)
        obs <- as.matrix(QUERY$Y[, c("psd1", "psd2", "psd3")])
        
        # Calcul R2
        r2_vec <- sapply(1:ncol(obs), function(j){
          cor(obs[,j], pred[,j], use = "complete.obs")^2
        })
        
        # Calcul RMSE
        rmse_vec <- sapply(1:ncol(obs), function(j){
          sqrt(mean((obs[,j] - pred[,j])^2, na.rm = TRUE))
        })
        
        MODEL[[i]][["eval"]][["R2"]] <- round(mean(r2_vec), 3)
        MODEL[[i]][["eval"]][["R2_sd"]] <- round(sd(r2_vec), 3)
        MODEL[[i]][["eval"]][["RMSE"]] <- round(mean(rmse_vec), 3)
        MODEL[[i]][["eval"]][["RMSE_sd"]] <- round(sd(rmse_vec), 3)
        
      } else if (i == "MLP") {
        # --- 1.2 SPECIFIQUE MLP (nnet) : R2 basée sur les prédictions Out-Of-Fold (OOF)
        # Utilisation des prédictions issues de la CV (générées dans model_continuous)
        # pour garantir une évaluation non biaisée.
        pred <- MODEL[[i]][["oof_preds"]]
        obs  <- MODEL[[i]][["oof_obs"]]
        
        # Calcul R2
        r2_vec <- sapply(1:ncol(obs), function(j){
          cor(obs[,j], pred[,j], use = "complete.obs")^2
        })
        
        # Calcul RMSE
        rmse_vec <- sapply(1:ncol(obs), function(j){
          sqrt(mean((obs[,j] - pred[,j])^2, na.rm = TRUE))
        })
        
        MODEL[[i]][["eval"]][["R2"]] <- round(mean(r2_vec), 3)
        MODEL[[i]][["eval"]][["R2_sd"]] <- round(sd(r2_vec), 3)
        MODEL[[i]][["eval"]][["RMSE"]] <- round(mean(rmse_vec), 3)
        MODEL[[i]][["eval"]][["RMSE_sd"]] <- round(sd(rmse_vec), 3)
        
      } else {
        # --- 1.3. Load final model data (Tidymodels)
        model_data <- lapply(1:(length(MODEL[[i]][["final_fit"]])), function(x){
          
          final_fit <- MODEL[[i]][["final_fit"]][[x]] %>%
            collect_predictions()
          
          # --- MULTI-OUTPUT : extraction obs/pred
          obs <- final_fit[, c("psd1","psd2","psd3")]
          
          pred_cols <- grep("^\\.pred", colnames(final_fit), value = TRUE)
          pred <- final_fit[, pred_cols]
          
          df <- data.frame(
            psd1 = obs[,1],
            psd2 = obs[,2],
            psd3 = obs[,3],
            pred1 = pred[,1],
            pred2 = pred[,2],
            pred3 = pred[,3]
          )
          
          return(df)
        }) %>% bind_rows()
        
        # --- Compute R-squared and RMSE into MODELS object
        obs <- as.matrix(model_data[, c("psd1","psd2","psd3")])
        pred <- as.matrix(model_data[, c("pred1","pred2","pred3")])
        
        # Calcul R2
        r2_vec <- sapply(1:ncol(obs), function(j){
          cor(obs[,j], pred[,j], use = "complete.obs")^2
        })
        
        # Calcul RMSE
        rmse_vec <- sapply(1:ncol(obs), function(j){
          sqrt(mean((obs[,j] - pred[,j])^2, na.rm = TRUE))
        })
        
        # Stockage
        MODEL[[i]][["eval"]][["R2"]] <- round(mean(r2_vec), 3)
        MODEL[[i]][["eval"]][["R2_sd"]] <- round(sd(r2_vec), 3)
        MODEL[[i]][["eval"]][["RMSE"]] <- round(mean(rmse_vec), 3)
        MODEL[[i]][["eval"]][["RMSE_sd"]] <- round(sd(rmse_vec), 3)
        # --- Memory cleanup after each iteration
        rm(model_data, df)
        gc()
      }
    } # for each model loop
  
  # --- 2. Variable importance - algorithm level
  # --- 2.1. Initialize function
  var_imp <- NULL
  
  # --- 2.1.2. General features - used later for plots
  features <- QUERY[["FOLDS"]][["resample_split"]][["splits"]][[1]]$data %>%
    dplyr::select(all_of(QUERY$SUBFOLDER_INFO$ENV_VAR))
  
  # --- 2.1.3. Loop over models
  for(i in MODEL$MODEL_LIST){
    
    if(i == "RF"){
      # --- 2.2 SPECIFIQUE RF : Extraction de la VIP native OOB
      rf_mod <- MODEL[[i]][["rfsrc_model"]]
      targets_lower <- c("psd1", "psd2", "psd3")
      
      out_list <- lapply(1:3, function(idx){
        imp_vals <- rf_mod$regrOutput[[ targets_lower[idx] ]]$importance
        
        # Sécurité : mettre les importances négatives (bruit) à 0
        imp_vals <- pmax(imp_vals, 0)
        
        data.frame(permutation = 1,
                   variable = names(imp_vals),
                   value = as.numeric(imp_vals),
                   cv = 1,
                   target = targets_lower[idx])
      })
      
      var_imp[[i]][["Raw"]] <- bind_rows(out_list)
      
    } else if (i == "MLP") {
      # --- 2.3 SPECIFIQUE MLP : VIP par permutation globale via DALEX
      mlp_mod <- MODEL[[i]][["nnet_model"]]
      targets_lower <- c("psd1", "psd2", "psd3")
      features_full <- QUERY$X[, QUERY$SUBFOLDER_INFO$ENV_VAR]
      
      out_list <- lapply(1:3, function(idx){
        target_name <- targets_lower[idx]
        
        # SÉCURITÉ 1 : S'assurer que y est un vecteur numérique pur (évite le rejet par DALEX)
        target_val <- as.numeric(unlist(QUERY$Y[, target_name]))
        
        # SÉCURITÉ 2 : Utiliser l'index [idx] au lieu du nom de colonne 
        # pour pallier à une potentielle perte de colnames() par predict.nnet()
        pred_fun <- function(m, d) {
          as.numeric(predict(m, newdata = d)[, idx])
        }
        
        explainer <- DALEX::explain(model = mlp_mod,
                                    data = features_full,
                                    y = target_val,
                                    predict_function = pred_fun,
                                    label = paste("MLP", target_name),
                                    verbose = FALSE)
        
        tmp <- DALEX::model_parts(explainer = explainer,
                                  loss_function = DALEX::loss_root_mean_square) %>%
          dplyr::filter(permutation != 0) %>%
          dplyr::filter(variable != "_baseline_") %>%
          group_by(permutation) %>%
          mutate(value = dropout_loss - dropout_loss[variable == "_full_model_"]) %>%
          dplyr::filter(variable != "_full_model_") %>% 
          ungroup() %>% 
          mutate(cv = 1,
                 target = target_name)
        
        return(tmp)
      })
      
      var_imp[[i]][["Raw"]] <- bind_rows(out_list)
      
    } else {
      # --- 2.4. Loop over the cross-validations (Tidymodels)
      var_imp[[i]][["Raw"]] <- lapply(1:CALL$NFOLD, function(x){
        
        id <- QUERY[["FOLDS"]][["resample_split"]][["splits"]][[x]][["in_id"]]
        
        features_x <- QUERY[["FOLDS"]][["resample_split"]][["splits"]][[x]]$data[id,] %>%
          dplyr::select(all_of(QUERY$SUBFOLDER_INFO$ENV_VAR))
        
        targets <- c("psd1","psd2","psd3")
        
        out_list <- lapply(targets, function(tgt){
          
          target <- QUERY[["FOLDS"]][["resample_split"]][["splits"]][[x]]$data[id,] %>%
            dplyr::select(all_of(tgt))
          
          m <- extract_fit_parsnip(MODEL[[i]][["final_fit"]][[x]])
          
          explainer <- explain_tidymodels(model = m,
                                          data = features_x,
                                          y = target)
          
          rm(target, m)
          gc()
          
          tmp <- model_parts(
            explainer = explainer,
            loss_function = loss_root_mean_square
          ) %>%
            dplyr::filter(permutation != 0) %>%
            dplyr::filter(variable != "_baseline_") %>%
            group_by(permutation) %>%
            mutate(value = dropout_loss - dropout_loss[variable == "_full_model_"]) %>%
            dplyr::filter(variable != "_full_model_") %>% 
            ungroup() %>% 
            mutate(cv = x,
                   target = tgt)
          
          return(tmp)
        })
        
        out <- bind_rows(out_list)
        return(out)
      }) %>% bind_rows() 
    }
    
    # --- 2.5. Further compute it as percentage for model-level plot
    if(sum(abs(var_imp[[i]][["Raw"]][["value"]]), na.rm = TRUE) > 0){
      var_imp[[i]][["Percent"]] <- var_imp[[i]][["Raw"]] %>%
        mutate(value = pmax(value, 0)) %>%
        group_by(cv, permutation, target) %>%
        mutate(value = if(sum(value) == 0) 0 else (value / sum(value) * 100)) %>%
        ungroup() %>%
        dplyr::select(variable, value, target) %>% 
        dplyr::filter(!is.na(value)) %>%
        mutate(variable = fct_reorder(variable, value, .desc = TRUE))
    } else {
      var_imp[[i]][["Percent"]] <- var_imp[[i]][["Raw"]]
    }
    
    # --- 2.6. Compute cumulative variable importance
    MODEL[[i]][["eval"]][["CUM_VIP"]] <- var_imp[[i]][["Percent"]] %>%
      group_by(target, variable) %>%
      summarise(average = mean(value), .groups = "drop") %>%
      group_by(target) %>%
      slice_max(order_by = average, n = 3) %>%
      summarise(total = sum(average)) %>%
      pull(total) %>%
      mean()
    
    # --- 2.7. Save row VIP
    MODEL[[i]][["vip"]] <- var_imp[[i]][["Percent"]]
    
  } # for each model loop
  
  # --- 3. Removing low quality algorithms
  for(i in MODEL$MODEL_LIST){
    if(MODEL[[i]][["eval"]][["R2"]] < 0.25 | is.na(MODEL[[i]][["eval"]][["R2"]])){
      MODEL$MODEL_LIST <- MODEL$MODEL_LIST[MODEL$MODEL_LIST != i]
      message(paste("--- EVAL : discarded", i, "due to R2 =", MODEL[[i]][["eval"]][["R2"]], "< 0.25 \n"))
    }
    
    if(MODEL[[i]][["eval"]][["CUM_VIP"]] < 50 | is.na(MODEL[[i]][["eval"]][["CUM_VIP"]])){
      MODEL$MODEL_LIST <- MODEL$MODEL_LIST[MODEL$MODEL_LIST != i]
      message(paste("--- EVAL : discarded", i, "due to CUM_VIP =", MODEL[[i]][["eval"]][["CUM_VIP"]], "< 50% \n"))
    }
  }
  
  # --- 4. Variable importance - Ensemble level
  if(CALL$ENSEMBLE == TRUE & length(MODEL$MODEL_LIST) > 1){
    ens_imp <- NULL
    message("--- VAR IMPORTANCE : compute ensemble")
    
    for(i in MODEL$MODEL_LIST){
      ens_imp <- rbind(ens_imp, var_imp[[i]][["Raw"]])
    } 
    
    # --- 4.2. Further compute it as percentage
    var_imp[["ENSEMBLE"]][["Percent"]] <- ens_imp %>%
      mutate(value = pmax(value, 0)) %>%
      group_by(cv, permutation, target) %>%
      mutate(value = if(sum(value) == 0) 0 else (value / sum(value) * 100)) %>%
      ungroup() %>%
      dplyr::select(variable, value, target) %>%
      mutate(variable = fct_reorder(variable, value, .desc = TRUE))
  } 
  
  # --- 5. Build the ensemble QC if there is an ensemble
  if(CALL$ENSEMBLE == TRUE & (length(MODEL$MODEL_LIST) > 1)){
    MODEL[["ENSEMBLE"]][["eval"]][["R2"]] <- lapply(MODEL$MODEL_LIST,
                                                    FUN = function(x){
                                                      x <- MODEL[[x]]$eval$R2
                                                    }) %>%
      unlist() %>% mean()
    
    MODEL[["ENSEMBLE"]][["eval"]][["CUM_VIP"]] <- var_imp[["ENSEMBLE"]][["Percent"]]  %>%
      group_by(variable) %>%
      summarise(average = mean(value), .groups = "drop") %>%
      dplyr::slice(1:3) %>%
      dplyr::select(average) %>%
      sum()
    
    MODEL[["ENSEMBLE"]][["vip"]] <- var_imp[["ENSEMBLE"]][["Percent"]]
  }
  
  # --- 6. Variable importance - Plot
  if(CALL$FAST == FALSE){
    plot_display <- CALL$HP$MODEL_LIST
  } else if(length(MODEL$MODEL_LIST) >= 1){
    plot_display <- MODEL$MODEL_LIST
  } else {
    plot_display <- NULL
  }
  
  if(!is.null(plot_display)){
    pdf(PDF_PATH, width = max(7, length(plot_display) * 3), height = 9)
    par(mfcol = c(3, length(plot_display)), mar = c(6,2,3,12))
    
    targets <- c("psd1", "psd2", "psd3")
    
    # --- 6.3.1. Algorithm level plot
    for(i in plot_display){
      if(i %in% MODEL$MODEL_LIST == TRUE){
        pal <- "#1F867B"
      } else {
        pal <- "#B64A60"
      }
      
      tmp <- var_imp[[i]][["Percent"]]
      
      for(tgt in targets){
        tmp_t <- tmp %>% dplyr::filter(target == tgt)
        
        # Eviter le crash si une target est vide pour une raison quelconque
        if(nrow(tmp_t) == 0) {
          plot.new()
          next
        }
        
        boxplot(tmp_t$value ~ tmp_t$variable, axes = FALSE, horizontal = TRUE,
                main = paste("VIP (", i, "-", tgt, ")"),
                sub = paste("R2 =", round(MODEL[[i]]$eval$R2, 2),
                            "| Cum. VIP(%) =", round(MODEL[[i]]$eval$CUM_VIP, 0)),
                col = pal, ylab = "", xlab = "Variable importance (%)")
        
        axis(side = 4, at = 1:nlevels(tmp_t$variable),
             labels = levels(tmp_t$variable), las = 2, cex.axis = 0.6)
        axis(side = 1, at = seq(0, 100, 10), labels = seq(0, 100, 10))
        abline(v = seq(0, 100, 10), lty = "dotted")
        box()
        box("figure", col="black", lwd = 1)
      }
    }
    
    # --- 6.3.2. Ensemble level plot
    if(CALL$ENSEMBLE == TRUE & (length(MODEL$MODEL_LIST) > 1)){
      tmp <- var_imp[["ENSEMBLE"]][["Percent"]]
      boxplot(tmp$value ~ tmp$variable, axes = FALSE, horizontal = TRUE,
              main = "PREDICTOR IMPORTANCE ( Ensemble )",
              sub = paste("Predictive performance (R2) =", round(MODEL[["ENSEMBLE"]]$eval$R2, 2), "; Cumulated var. importance (%; top 3) =", round(MODEL[["ENSEMBLE"]]$eval$CUM_VIP, 0)),
              col = "gray50", ylab = "", xlab = "Variable importance (%)")
      
      # Utilisation de nlevels() au lieu de ncol() pour correspondre exactement à l'axe
      axis(side = 4, at = 1:nlevels(tmp$variable), labels = levels(tmp$variable), las = 2, cex.axis = 0.6)
      axis(side = 1, at = seq(0, 100, 5), labels = seq(0, 100, 5))
      abline(v = seq(0, 100, 10), lty = "dotted")
      box()
      box("figure", col="black", lwd = 1)
    }
    
    dev.off()
  }
  
  return(MODEL)
  
} # END FUNCTION
