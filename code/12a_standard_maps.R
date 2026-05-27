#' =============================================================================
#' @name standard_maps
#' @description Simple function for building standard mean and SD maps from
#' projection data
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @param MONTH list corresponding to the groups of month to computes maps from.
#' @return plots mean and uncertainty maps per model or ensemble

standard_maps <- function(FOLDER_NAME = NULL, SUBFOLDER_NAME = NULL,
                          MONTH = list(c(10,11,12,1,2,3), 4:9)){

  # --- 1. Initialize function
  set.seed(123)  # Set seed for reproducibility
  
  # --- 1.1. Start logs - append file
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/log.txt"), open = "a"),
                       START = TRUE)
  message(paste(Sys.time(), "******************** START : standard_maps ********************"))
  
  # --- 1.2. Parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData"))
  
  # --- 1.3. Check for projections
  if((length(MODEL$MODEL_LIST) == 0) & CALL$FAST == TRUE){
    message("No validated algorithms to display projections from")
    log_sink(FILE = sinkfile, START = FALSE)
    return(NULL)
  }
  
  # --- 1.3. Initialize loop over algorithm or targets for later
  if(CALL$DATA_TYPE == "proportions"){loop_over <- 1:ncol(QUERY$Y)
  }else if(CALL$FAST == FALSE){loop_over <- CALL$HP$MODEL_LIST
  }else {loop_over <- MODEL$MODEL_LIST}
  
  # --- 1.4. Create PDF saving (CORRECTION : fermeture sécurisée via on.exit)
  pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/",SUBFOLDER_NAME,"/05_standard_maps.pdf"))
  on.exit({
    if(names(dev.cur()) != "null device") dev.off()
  }, add = TRUE)
  
  # --- 1.5. Set initial plot layout & requirements
  par(mfrow = c(4,3), mar = c(2,2,7,1))
  
  # --- 1.6. Set base raster
  CALL$ENV_DATA <- lapply(CALL$ENV_DATA, function(x) terra::rast(x)) # Unpack the rasters first
  QUERY$MESS <- terra::rast(QUERY$MESS) # Unpack the mess as well
  r0 <- CALL$ENV_DATA[[1]][[1]]
  
  # --- 1.7 Land mask
  land <- r0
  land[is.na(land)] <- 9999
  land[land != 9999] <- NA
  
  # --- 2. Plot the quality checks
  # --- 2.1. Compute the recommendation table
  MODEL$recommandations <- qc_recommandations(QUERY = QUERY, MODEL = MODEL, DATA_TYPE = CALL$DATA_TYPE)
  
  traffic_col <- rep(MODEL$recommandations$COL, each = 4)
  traffic_val <- MODEL$recommandations[,1:4] %>% as.matrix() %>% t() %>% c()
  traffic_col[which(as.numeric(traffic_val) == 0)] <- "white"
  
  # --- 2.2. Plot the traffic lights and recommendations
  plot.new()
  if(CALL$DATA_TYPE != "proportions"){
    mtext(paste("QUALITY CHECK \n",
                QUERY$annotations$scientificname, 
                "\n ID:", QUERY$annotations$worms_id))
  } else {
    mtext("QUALITY CHECK \n for proportions")
  }
  
  par(mar = c(1,3,7,1), xpd = NA)
  plot(x = rep(1:4, nrow(MODEL$recommandations)), 
       y = rep(nrow(MODEL$recommandations):1, each = 4), 
       axes = FALSE, cex = 3,
       xlim = c(0,5), ylim = c(0,nrow(MODEL$recommandations)+1), 
       ylab = "", xlab = "",
       pch = 21, col = "black", 
       bg = traffic_col)
  axis(side = 3, at = 1:4, labels = c("A priori \n var. imp.","Perdictive \n performance","Cumulative \n var. imp.","Projection \n uncertainty"), tick = FALSE, line = NA, cex.axis = 1, las = 2)
  axis(side = 2, at = nrow(MODEL$recommandations):1, labels = rownames(MODEL$recommandations), tick = FALSE, line = NA, las = 2, cex.axis = 0.7)
  axis(side = 4, at = nrow(MODEL$recommandations):1, labels = MODEL$recommandations$Recommandation, tick = FALSE, line = NA, las = 2, cex.axis = 0.7)
  abline(h = 0) # PDF wide separator
  plot.new()
  par(mar = c(2,2,4,1))
  
  # --- 3. Plot the legends
  hsi_pal <- inferno_pal(100)
  plot.new()
  
  # --- 3.1.1. Extract the plot true scale (Kept for SD / Uncertainty maps)
  if(CALL$DATA_TYPE == "continuous"){
    plot_scale <- lapply(loop_over, FUN = function(z){
      lapply(MODEL[[z]]$proj$y_hat, function(tgt) {
        apply(tgt, 1, function(x) mean(x, na.rm = TRUE))
      }) %>% unlist()
    }) %>% unlist() %>% quantile(0.95, na.rm = TRUE)
    
    if(is.na(plot_scale) || plot_scale <= 0) plot_scale <- 1 
    
  } else {
    plot_scale <- 1
  }
  
  # --- 3.1.2. Plot the global dummy colorbar
  colorbar.plot(x = 0.5, y = 0, strip = seq(0,1,length.out = 100),
                strip.width = 0.3, strip.length = 2.7, col = hsi_pal, border = "black")
  axis(side = 1, at = seq(0, 1, length.out = 5), labels = round(seq(0, 1 * plot_scale, length.out = 5), 2))
  text(x = 0.5, y = 0.3, "Habitat Suitability Index", adj = 0.5)
  
  # --- 3.2. Observation vs 75% quartile
  plot.new()
  points(x = 0.1, y = 0.4, pch = 22, col = "black", bg = "gray80", cex = 5)
  text(x = 0.2, y = 0.4, "Q75 Habitat Suitability Index", pos = 4)
  points(x = 0.1, y = 0.1, pch = 20, col = "black", cex = 2)
  text(x = 0.2, y = 0.1, "Observation", pos = 4)
  
  # --- 3.3. MESS x SD 2D color scale
  par(mar = c(5,5,3,2), xpd = FALSE)
  bivar_pal <- colmat(nbreaks = 100)
  colmat_plot(bivar_pal, xlab = "Standard deviation", ylab = "MESS value")
  axis(side = 1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = round(seq(0, 1 * plot_scale*0.5, length.out = 6), 2))
  axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, -20, -40, -60, -80, -100), las = 2)
  par(mar = c(4,2.5,3,1))
  
  # --- 4. Build model-level outputs
  for(i in loop_over){
    if(CALL$DATA_TYPE == "proportions"){
      val_raw_list <- list(MODEL[["MBTR"]][["proj"]][["y_hat"]][,,i,]) 
      model_name <- "MBTR"
    } else {
      val_raw_list <- MODEL[[i]][["proj"]][["y_hat"]]
      model_name <- i
    }
    
    for(t in seq_along(val_raw_list)){
      val_raw <- val_raw_list[[t]]
      tgt_name <- if(CALL$DATA_TYPE == "proportions") colnames(QUERY$Y)[i] else colnames(QUERY$Y)[t]
      
      if(CALL$DATA_TYPE == "proportions"){
        obs_vec <- QUERY$Y[[i]] 
      } else {
        obs_vec <- QUERY$Y[[t]] 
      }
      
      for(j in seq_along(MONTH)){
        
        # Mean value - CORRECTION : On ne divise plus par plot_scale pour avoir la vraie valeur de la cible
        val <- apply(val_raw[,,j], 1, function(x)(x = mean(x, na.rm = TRUE)))
        r_m <- terra::rast(r0, vals = val)
        
        # --- CORRECTION DYNAMIQUE : On calcule val_min et val_max APRES avoir r_m
        val_min <- min(c(min(obs_vec, na.rm = TRUE), min(terra::values(r_m, mat=FALSE), na.rm = TRUE))) 
        val_max <- max(c(max(obs_vec, na.rm = TRUE), max(terra::values(r_m, mat=FALSE), na.rm = TRUE)))
        
        # Coefficient of variation
        if(length(MONTH[[j]]) > 1){
          val <- apply(val_raw[,,MONTH[[j]]], c(1,3), function(x)(x = sd(x, na.rm = TRUE))) %>%
            apply(1, function(x)(x = mean(x, na.rm = TRUE)))
        } else {
          val <- apply(val_raw[,,j], 1, function(x)(x = sd(x, na.rm = TRUE)))
        }
        
        r_sd <- terra::rast(r0, vals = val)
        r_sd[r_sd > plot_scale*25] <- plot_scale*25
        r_sd[r_sd <= 0] <- 1e-10 
        
        # MESS
        r_mess <- QUERY$MESS[[j]]*-1
        if(nlyr(r_mess) > 1){r_mess <- app(r_mess, mean, na.rm = TRUE)}
        r_mess[r_mess<0] <- 1e-10 
        r_mess[r_mess>100] <- 100 
        
        # --- 4.3. Plot maps (CORRECTION: legend = TRUE et ajout de l'argument 'range')
        if(CALL$DATA_TYPE == "proportions"){
          tmp_idx <- which(QUERY$annotations$worms_id == tgt_name)
          plot(r_m, col = hsi_pal, range = c(val_min, val_max), legend=TRUE, cex.main = 1,
               main = paste("Projection for", QUERY$annotations$scientificname[tmp_idx], "\n Month:", paste(MONTH[[j]], collapse = ",")))
          mtext(text = paste("Predictive performance (", names(MODEL[["MBTR"]][["eval"]])[1], ") =", MODEL[["MBTR"]][["eval"]][[1]],
                             "\n Colorbar scale:", format(round(terra::minmax(r_m)[2], 5), scientific = TRUE)),
                side = 1, line = 3, cex = 0.6)
        } else {
          title_text <- if(length(val_raw_list) > 1) paste("Projection (", model_name, "-", tgt_name, ") \n Month:", paste(MONTH[[j]], collapse = ",")) 
          else paste("Projection (", model_name, ") \n Month:", paste(MONTH[[j]], collapse = ","))
          
          # La bidouille [m_min:m_max] est supprimée au profit de la vraie gestion avec 'range'
          plot(r_m, col = hsi_pal, range = c(val_min, val_max), legend=TRUE, cex.main = 1, main = title_text)
          mtext(text = paste("Predictive performance (", names(MODEL[[model_name]][["eval"]])[1], ") =", MODEL[[model_name]][["eval"]][[1]]),
                side = 1, line = 2, cex = 0.7)
        }
        
        plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)
        box("figure", col="black", lwd = 1)
        
        # Observations
        plot(r_m > quantile(terra::values(r_m, mat=FALSE), 0.75, na.rm = TRUE), col = c("white","gray80"), legend = FALSE, main = "Observations")
        plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)
        box("figure", col="black", lwd = 1)
        
        tmp_obs <- QUERY$S
        obs_vec_sub <- obs_vec
        
        # CORRECTION DYNAMIQUE : Les couleurs des points s'alignent sur val_min et val_max
        if(CALL$DATA_TYPE == "continuous"){
          points(tmp_obs$decimallongitude, tmp_obs$decimallatitude,
                 col = col_numeric("inferno", domain = c(val_min, val_max), alpha = 0.2, na.color = hsi_pal[100])(obs_vec_sub), pch = 20, cex = 0.6)
        } else if(CALL$DATA_TYPE == "proportions") {
          points(tmp_obs$decimallongitude, tmp_obs$decimallatitude,
                 col = col_numeric("inferno", domain = c(val_min, val_max), alpha = 0.2)(obs_vec_sub), pch = 20, cex = 0.6)
        } else {
          points(tmp_obs$decimallongitude, tmp_obs$decimallatitude, col = "black", pch = 20, cex = 0.6)
        }
        
        # Uncertainties
        plot(land, col = "antiquewhite4", legend=FALSE, main = "Uncertainties")
        box("figure", col="black", lwd = 1)
        
        max_sd <- max(terra::values(r_sd, mat=FALSE), na.rm = TRUE)
        if(is.infinite(max_sd) || is.na(max_sd) || max_sd <= 0) max_sd <- 1e-5
        
        r <- bivar_map(rasterx = r_sd, rastery = r_mess, colormatrix = bivar_pal,
                       cutx = seq(0, max_sd, length.out = 101), cuty = 0:100)
        plot(r[[1]], col = r[[2]], legend=FALSE, add = TRUE)
        
        eval_mod <- if(CALL$DATA_TYPE == "proportions") "MBTR" else model_name
        mtext(text = paste("Projection uncertainty (", names(MODEL[[eval_mod]][["eval"]])[3], ") =", round(MODEL[[eval_mod]][["eval"]][[3]],2)),
              side = 1, line = 2, cex = 0.7)
        
      } # End j month loop
    } # End t target loop
  } # End i model loop
  
  # --- 5. Build ensemble quality checks
  if(CALL$ENSEMBLE == TRUE & (length(MODEL$MODEL_LIST) > 1)){
    rec <- qc_recommandations(QUERY = QUERY, MODEL = MODEL, DATA_TYPE = CALL$DATA_TYPE, ENSEMBLE = TRUE)
    
    traffic_col <- rep(rec$COL, each = 4)
    traffic_val <- rec[,1:4] %>% as.matrix() %>% t() %>% c()
    traffic_col[which(as.numeric(traffic_val) == 0)] <- "white"
    
    plot.new()
    par(mar = c(1,3,7,1), xpd = NA)
    plot(x = rep(1:4, nrow(rec)), y = rep(nrow(rec):1, each = 4), axes = FALSE, cex = 3,
         xlim = c(0,5), ylim = c(0,nrow(rec)+1), ylab = "", xlab = "",
         pch = 21, col = "black", bg = traffic_col)
    axis(side = 3, at = 1:4, labels = c("A priori \n var. imp.","Perdictive \n performance","Cumulative \n var. imp.","Projection \n uncertainty"), tick = FALSE, line = NA, cex.axis = 1, las = 2)
    axis(side = 2, at = nrow(rec):1, labels = rownames(rec), tick = FALSE, line = NA, las = 2, cex.axis = 0.7)
    axis(side = 4, at = nrow(rec):1, labels = rec$Recommandation, tick = FALSE, line = NA, las = 2, cex.axis = 0.7)
    plot.new()
    par(mar = c(4,2.5,3,1))
  } # End ensemble QC
  
  # --- 6. Build ensemble outputs
  if(CALL$ENSEMBLE == TRUE & (length(MODEL$MODEL_LIST) > 1)){
    
    n_targets <- length(MODEL[[ MODEL$MODEL_LIST[1] ]][["proj"]][["y_hat"]])
    y_ens_list <- list()
    
    for(t in 1:n_targets){
      tgt_ens <- NULL
      for(i in MODEL$MODEL_LIST){
        tgt_ens <- abind::abind(tgt_ens, MODEL[[i]][["proj"]][["y_hat"]][[t]], along = 2)
      }
      y_ens_list[[t]] <- tgt_ens
    }
    
    r_mess_ens <- app(QUERY$MESS * -1, mean, na.rm = TRUE)
    r_mess_ens[r_mess_ens < 0] <- 1e-10
    r_mess_ens[r_mess_ens > 100] <- 100
    
    for(t in seq_along(y_ens_list)){
      val_raw <- y_ens_list[[t]]
      tgt_name <- colnames(QUERY$Y)[t]
      
      # CORRECTION DYNAMIQUE ENSEMBLE : valeurs pures, limites calculées après
      val <- apply(val_raw, 1, function(x)(x = mean(x, na.rm = TRUE)))
      r_m <- terra::rast(r0, vals = val) 
      
      obs_vec <- QUERY$Y[[t]]
      tmp_obs <- QUERY$S[which(obs_vec > 0), ]
      obs_vec_sub <- obs_vec[which(obs_vec > 0)]
      
      val_min <- min(c(min(obs_vec, na.rm = TRUE), min(terra::values(r_m, mat=FALSE), na.rm = TRUE)))
      val_max <- max(c(max(obs_vec, na.rm = TRUE), max(terra::values(r_m, mat=FALSE), na.rm = TRUE)))
      
      val_sd <- apply(val_raw, 1, function(x)(x = sd(x, na.rm = TRUE)))
      r_sd <- terra::rast(r0, vals = val_sd)
      r_sd[r_sd > plot_scale*0.25] <- plot_scale*0.25
      r_sd[r_sd <= 0] <- 1e-10 
      
      title_text <- if(n_targets > 1) paste("Projection ( Ensemble -", tgt_name, ")") else "Projection ( Ensemble )"
      
      plot(r_m, col = hsi_pal, range = c(val_min, val_max), legend=TRUE, main = title_text)
      mtext(text = paste("Predictive performance (", names(MODEL[["ENSEMBLE"]][["eval"]])[1], ") =", round(MODEL[["ENSEMBLE"]][["eval"]][[1]],2)),
            side = 1, line = 2, cex = 0.7)
      plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)
      box("figure", col="black", lwd = 1)
      
      plot(r_m > quantile(terra::values(r_m, mat=FALSE), 0.75, na.rm = TRUE), col = c("white","gray80"), legend=FALSE, main = "Observations")
      plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)
      box("figure", col="black", lwd = 1)
      
      if(CALL$DATA_TYPE == "continuous"){
        points(tmp_obs$decimallongitude, tmp_obs$decimallatitude,
               col = col_numeric("inferno", domain = c(val_min, val_max), alpha = 0.2, na.color = hsi_pal[100])(obs_vec_sub), pch = 20, cex = 0.6)
      } else {
        points(tmp_obs$decimallongitude, tmp_obs$decimallatitude, col = "black", pch = 20)
      }
      
      plot(land, col = "antiquewhite4", legend=FALSE, main = "Uncertainties")
      box("figure", col="black", lwd = 1)
      
      max_sd_ens <- max(terra::values(r_sd, mat=FALSE), na.rm = TRUE)
      if(is.infinite(max_sd_ens) || is.na(max_sd_ens) || max_sd_ens <= 0) max_sd_ens <- 1e-5
      
      r <- bivar_map(rasterx = r_sd, rastery = r_mess_ens, colormatrix = bivar_pal,
                     cutx = seq(0, max_sd_ens, length.out = 101), cuty = 0:100)
      plot(r[[1]], col = r[[2]], legend=FALSE, add = TRUE)
      mtext(text = paste("Projection uncertainty (", names(MODEL[["ENSEMBLE"]][["eval"]])[3], ") =", round(MODEL[["ENSEMBLE"]][["eval"]][[3]],2)),
            side = 1, line = 2, cex = 0.7)
      
    } # End t target loop for Ensemble
  } # End if ENSEMBLE = TRUE
  
  # --- 7. Wrap up and save
  log_sink(FILE = sinkfile, START = FALSE)
  save(MODEL, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData"),
       compress = "gzip", compression_level = 6)
  return(SUBFOLDER_NAME)
  
} # END FUNCTION
