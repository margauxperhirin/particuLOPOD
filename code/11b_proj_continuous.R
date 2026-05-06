#' =============================================================================
#' @name proj_continuous
#' @description computes spatial projections for the continuous data sub-pipeline
#' @param CALL the call object from the master pipeline
#' @param QUERY the query object from the master pipeline
#' @param MODEL the models object from the master pipeline
#' @return an updated model list object containing the projections objects
#' embedded in each model sub-list.

proj_continuous <- function(QUERY, MODEL, CALL){

  # --- 1. Initialize function
  # --- 1.1. Open base raster and values
  CALL$ENV_DATA <- lapply(CALL$ENV_DATA, function(x) terra::rast(x)) # Unpack the rasters first
  r0 <- CALL$ENV_DATA[[1]][[1]] # Base raster 
  
  # --- 1.2. Define the projections to compute
  # Only the algorithm that passed QC if the argument FAST equals TRUE
  loop_over <- if(CALL$FAST) MODEL$MODEL_LIST else CALL$HP$MODEL_LIST
  
  # --- 2. Define bootstraps
  # --- 2.1. Target transformation
  if(CALL$DATA_TYPE == "continuous" & !is.null(CALL$TARGET_TRANSFORMATION)){
    message("--- PROJ : Transforming the target variable according to the provided function")
    source(CALL$TARGET_TRANSFORMATION) # Load the transformation function
    tmp <- target_transformation(QUERY$Y$measurementvalue, REVERSE = FALSE) # Transformation
    Y <- data.frame(target_transformation(QUERY$Y$measurementvalue, REVERSE = FALSE)$out) # New target
    colnames(Y) <- "measurementvalue"
    QUERY[["target_transformation"]][["yj_obj"]] <- tmp$yj_obj # Save the transformation parameters
  } else {
    Y <- QUERY$Y
  }
  
  # --- 2.2. Re-assemble all query tables
  tmp <- cbind(Y, QUERY$X)
  
  # --- 2.3. Run the bootstrap generation from tidy models
  boot_split <- bootstraps(tmp, times = CALL$N_BOOTSTRAP)
  rm(tmp) # Intermediate cleanup (keep Y to extract some names)
  gc()
  
  # --- 3. Start the loop over algorithms
  for(i in loop_over){
    # --- 3.1. Fit model on bootstrap
    if(i == "RF"){
      # --- SPECIFIQUE RF : Entraînement manuel sur chaque bootstrap
      message("--- PROJ : Fitting RF on bootstraps...")
      
      target_names <- colnames(Y)
      # Détection auto : régression simple ou multivariée
      if(length(target_names) > 1) {
        frm <- as.formula(paste0("Multivar(", paste(target_names, collapse = ", "), ") ~ ."))
        is_multivar <- TRUE
      } else {
        frm <- as.formula(paste0(target_names, " ~ ."))
        is_multivar <- FALSE
      }
      
      boot_fit_rf <- lapply(boot_split$splits, function(split) {
        boot_data <- rsample::analysis(split) %>%
          dplyr::select(all_of(c(target_names, QUERY$SUBFOLDER_INFO$ENV_VAR)))
        
        # Entraînement sans importance des variables pour gagner du temps en projection
        randomForestSRC::rfsrc(frm, data = as.data.frame(boot_data), ntree = 500, importance = "none")
      })
      
    } else if (i == "MLP") {
      # --- SPECIFIQUE MLP : Entraînement manuel sur chaque bootstrap
      message("--- PROJ : Fitting MLP (nnet) on bootstraps...")
      require(nnet)
      
      target_names <- colnames(Y)
      mlp_formula <- as.formula(paste0("cbind(", paste(target_names, collapse = ", "), ") ~ ."))
      best_params <- MODEL[[i]][["best_params"]]
      
      boot_fit_mlp <- lapply(boot_split$splits, function(split) {
        boot_data <- rsample::analysis(split) %>%
          dplyr::select(all_of(c(target_names, QUERY$SUBFOLDER_INFO$ENV_VAR)))
        
        for(v in QUERY$SUBFOLDER_INFO$ENV_VAR) {
          boot_data[[v]] <- (boot_data[[v]] - QUERY$X_norm_params$means[v]) / QUERY$X_norm_params$sds[v]
        }
        
        nnet::nnet(formula = mlp_formula,
                   data = boot_data,
                   size = best_params$hidden_units,
                   epochs = params$epochs,
                   linout = TRUE,
                   trace = FALSE,
                   maxit = 1000)
      })
      
    } else {
      # --- TIDYMODELS : fit_resamples() does not save models by default. Thus the control_resamples()
      boot_fit <- MODEL[[i]][["final_wf"]] %>%
        fit_resamples(resamples = boot_split,
                      control = control_resamples(extract = function (x) extract_fit_parsnip(x),
                                                  verbose = FALSE)) %>%
        unnest(.extracts)
    }
    
    # --- 4. Loop over month for predictions
    y_hat <- NULL
    for(m in seq_along(CALL$ENV_DATA)){
      
      # --- 4.1. Load the right features
      features <- terra::subset(CALL$ENV_DATA[[m]], QUERY$SUBFOLDER_INFO$ENV_VAR) %>% 
        terra::as.data.frame()
      
      # --- 4.2. Compute one prediction per bootstrap
      if(i == "RF"){
        # --- SPECIFIQUE RF : Prédiction et formatage mimant tidymodels
        boot_proj_list <- lapply(boot_fit_rf, function(mod) {
          p <- predict(mod, newdata = as.data.frame(features))
          
          if (is_multivar) {
            pred_df <- lapply(target_names, function(tgt) {
              p$regrOutput[[tgt]]$predicted
            }) %>% as.data.frame()
            colnames(pred_df) <- paste0(".pred_", target_names)
          } else {
            pred_df <- data.frame(.pred = p$predicted)
            colnames(pred_df) <- paste0(".pred_", target_names)
          }
          return(pred_df)
        })
        
        # On recrée l'objet attendu par la suite du code
        boot_proj <- tibble(
          id = boot_split$id,
          proj = boot_proj_list
        )
        
      } else if (i == "MLP") {
        # --- SPECIFIQUE MLP : Prédiction et formatage mimant tidymodels
        target_names <- colnames(Y)
        
        features_norm <- as.data.frame(features)
        
        for(v in QUERY$SUBFOLDER_INFO$ENV_VAR) {
          features_norm[[v]] <- (features_norm[[v]] - QUERY$X_norm_params$means[v]) / QUERY$X_norm_params$sds[v]
        } # Normalisation before projection 
        
        boot_proj_list <- lapply(boot_fit_mlp, function(mod) {
          # predict.nnet sur de la multivariée renvoie une matrice
          p <- predict(mod, newdata = as.data.frame(features))
          pred_df <- as.data.frame(p)
          
          # Renommage pour coller au standard tidymodels attendu par l'étape 4.3
          colnames(pred_df) <- paste0(".pred_", target_names)
          return(pred_df)
        })
        
        # On recrée l'objet attendu
        boot_proj <- tibble(
          id = boot_split$id,
          proj = boot_proj_list
        )
        
      } else {
        # --- TIDYMODELS : As we extracted the model information in a supplementary column...
        boot_proj <- boot_fit %>%
          mutate(proj = purrr::map(.extracts, function(x){
            p <- predict(x, features)
            
            # récupérer toutes les colonnes .pred*
            pred_cols <- grep("^\\.pred", colnames(p), value = TRUE)
            return(p[, pred_cols, drop = FALSE])
          }))
      }
      
      # --- 4.3. First transform the object into a cell x bootstrap matrix
      # /!\ Need to create a unique row identifier for pivot_wider to work...
      boot_proj2 <- boot_proj %>%
        dplyr::select(id, proj) %>%
        unnest(c(id, proj)) %>%
        as.data.frame()
      
      # récupérer colonnes de prédiction
      pred_cols <- grep("^\\.pred", colnames(boot_proj2), value = TRUE)
      
      # créer une matrice (cellules x bootstrap x variables)
      boot_proj <- lapply(pred_cols, function(pc){
        tmp <- boot_proj2 %>%
          dplyr::select(id, all_of(pc)) %>%
          group_by(id) %>%
          mutate(row = row_number()) %>%
          pivot_wider(names_from = id, values_from = all_of(pc)) %>%
          dplyr::select(-row)
        
        return(as.matrix(tmp))
      })
      
      # --- 4.4. Assign the desired values to the non-NA cells in the base raster
      boot_proj <- lapply(boot_proj, function(mat){
        apply(mat, 2, function(x){
          r <- terra::values(r0)
          r[!is.na(r)] <- x
          return(r)
        })
      })
      
      # --- 4.5. Concatenate with previous month
      y_hat <- lapply(seq_along(boot_proj), function(k){
        abind(y_hat[[k]], boot_proj[[k]], along = 3)
      })
      
      message(paste("--- PROJ :", i,"- month", m, "done \t"))
      
      rm(boot_proj, features)
      gc()
    } # for m month
    
    # --- 4.6. Reverse transformation
    if(CALL$DATA_TYPE == "continuous" & !is.null(CALL$TARGET_TRANSFORMATION)){
      message("--- PROJ : reverse transformation of the target")
      y_hat <- apply(y_hat, -1, function(x){
        x <- target_transformation(x, REVERSE = TRUE, PARAM = QUERY$target_transformation)
      })
    } # end if transformation
    
    # --- 5. Cut spatial discontinuities
    if (!is.null(CALL$CUT)) {
      y_hat <- apply(y_hat, c(2, 3), function(x) {
        xy <- QUERY$S %>% dplyr::select(decimallongitude, decimallatitude) %>%
          .[QUERY$Y == 1, ] # Extract values
        x2 <- x
        x[x < (CALL$CUT * max(x, na.rm = TRUE))] <- 0 # Set threshold
        r_patch <- terra::patches(setValues(r0, x), zeroAsNA = TRUE) # Compute patches
        id_patch <- terra::extract(r_patch, xy)$patches %>% unique() %>% .[!is.na(.)] # Verify if each patch has observations
        
        if(length(id_patch) > 0) { # Make sure there is patches 
          r_patch[!(r_patch %in% id_patch)] <- 0 # Remove the ones that dont
          r_patch[r_patch > 0] <- 1 # scale
          x <- x2 * terra::values(r_patch)
        } else {
          x <- x2
          message(paste("--- CUT : could not been applied to this projection of", i, 
                        "because it has a maximum suitability index inferior to the cutting threshold \n"))
        }
        
        return(x)
      })
    } # if cut
    
    # --- 6. Compute the average CV across bootstrap runs as a QC
    if(dim(y_hat[[1]])[[2]] == CALL$N_BOOTSTRAP) { # Modification mineure : y_hat est une liste de matrices, vérif sur y_hat[[1]]
      NSD_vec <- sapply(y_hat, function(arr){
        nsd <- apply(arr, c(1,3), function(x) sd(x, na.rm = TRUE)) %>%
          mean(na.rm = TRUE)
        
        nsd / mean(arr, na.rm = TRUE)
      })
      
      NSD <- mean(NSD_vec)
      # NSD <- NSD/mean(y_hat, na.rm = TRUE) # i think i hqve to remove this also
    } else {
      NSD <- NA
      message(paste("--- PROJ: Model", i, " discarded, bootstrap did not complete"))
    }
    
    # --- 7. Append the MODEL object
    # --- 7.1. Save the evaluation metric and projections
    MODEL[[i]][["proj"]][["y_hat"]] <- y_hat
    MODEL[[i]][["eval"]][["NSD"]] <- NSD
    
    # --- 7.2. Discard low quality models according to NSD
    # Fixed at 50% average across the projection
    # if(MODEL[[i]][["eval"]][["NSD"]] > 0.5 | is.na(MODEL[[i]][["eval"]][["NSD"]])){
    #   MODEL$MODEL_LIST <- MODEL$MODEL_LIST[MODEL$MODEL_LIST != i]
    #   message(paste("--- EVAL : discarded", i, "due to NSD =", MODEL[[i]][["eval"]][["NSD"]], "> 0.5 \n"))
    # }
    
    # --- 7.3. Discard if the PRE_VIP QC in the query is too low
    # In case of FAST = TRUE; the species was kept until here but we do not want
    # to compute an ensemble for it.
    if(QUERY$eval$PRE_VIP < 0.05){
      MODEL$MODEL_LIST <- MODEL$MODEL_LIST[MODEL$MODEL_LIST != i]
      message(paste("--- EVAL : discarded", i, "due to PRE_VIP =", QUERY$eval$PRE_VIP, "< 0.05 \n"))
    }
    
  } # for i model loop
  
  # --- 8. Build the ensemble NSD if there is an ensemble
  if(CALL$ENSEMBLE == TRUE & (length(MODEL$MODEL_LIST) > 1)){
    MODEL[["ENSEMBLE"]][["eval"]][["NSD"]] <- lapply(MODEL$MODEL_LIST,
                                                     FUN = function(x){
                                                       x <- MODEL[[x]]$eval$NSD
                                                     }) %>%
      unlist() %>% mean()
  } # End ENSEMBLE QC
  
  return(MODEL)
  
} # END FUNCTION
