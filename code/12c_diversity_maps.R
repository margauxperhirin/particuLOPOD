#' =============================================================================
#' @name diversity_maps
#' @description Function extracting all projections by species and bootstrap.
#' Then computes a set of diversity metrics and plots
#' @param FOLDER_NAME name of the corresponding folder
#' @param HILL a value of Hill number to use for diversity computation
#' @param MONTH a list of month to concatenate together
#' @return plots mean and uncertainty maps per diversity metric
#' @return a diversity and mess object per projection, species and bootstrap.
#' @return a .nc matching the EMODnet standards with the same informations as 
#' above. - Saved in FOLDERNAME.

diversity_maps <- function(FOLDER_NAME = NULL,
                           HILL = c(0,1,2),
                           MONTH = list(c(10,11,12,1,2,3),
                                        4:9)){

  # --- 1. Global parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  
  # Baseline raster
  CALL$ENV_DATA <- lapply(CALL$ENV_DATA, function(x) terra::rast(x)) # Unpack the rasters first
  r0 <- CALL$ENV_DATA[[1]][[1]]
  
  # Land mask
  land <- r0
  land[is.na(land)] <- 9999
  land[land != 9999] <- NA
  
  # --- 2. Predict first - assemble later estimation
  message(paste0(Sys.time(), "--- DIVERSITY : build the ensembles - START"))
  
  # --- 2.1. Which subfolder list and which model in the ensemble
  all_files <- list.files(paste0(project_wd, "/output/", FOLDER_NAME), recursive = TRUE)
  model_files <- unique(dirname(all_files[grepl("MODEL.RData", all_files)])) %>% .[1:30]
  ensemble_files <- mclapply(model_files, function(x){
    memory_cleanup() # low memory use
    
    load(paste0(project_wd, "/output/", FOLDER_NAME,"/", x, "/MODEL.RData"))
    if(length(MODEL$MODEL_LIST) >= 1){
      load(paste0(project_wd, "/output/", FOLDER_NAME,"/", x, "/QUERY.RData"))
      return(list(SUBFOLDER_NAME = x, MODEL_LIST = MODEL$MODEL_LIST, Y = QUERY$Y, MESS = QUERY$MESS, REC = MODEL$recommandations))
    } else {
      return(NULL)
    } # if model list
  }, mc.cores = round(MAX_CLUSTERS/2, 0), mc.preschedule = FALSE) %>% 
    .[lengths(.) != 0] %>% 
    .[grep("Error", ., invert = TRUE)] # to exclude any API error or else
  
  # --- 2.2. Loop over the files
  message(paste0(Sys.time(), "--- DIVERSITY: build the ensembles - loop over files"))
  tmp <- mclapply(ensemble_files, function(x){
    memory_cleanup() # low memory use
    
    # --- 2.2.1. Load MODEL files
    load(paste0(project_wd, "/output/", FOLDER_NAME,"/", x$SUBFOLDER_NAME, "/MODEL.RData"))
    
    # --- 2.2.2. Extract projections in a matrix
    # If there is more than 1 algorithm, we extract and average across algorithm
    # Output matrix is cell x bootstrap x month
    if(length(x$MODEL_LIST) > 1){
      m <- lapply(x$MODEL_LIST, function(y){
        MODEL[[y]][["proj"]]$y_hat
      }) %>% abind(along = 4) %>% apply(c(1,2,3), function(z)(z = mean(z, na.rm = TRUE)))
    } else {
      m <- MODEL[[x$MODEL_LIST]][["proj"]]$y_hat
    } # end if
  }, mc.cores = round(MAX_CLUSTERS/2, 0), mc.cleanup = TRUE)
  
  # --- 2.3. Re-arrange the projections in an array
  message(paste0(Sys.time(), "--- DIVERSITY: build the ensembles - format to array"))
  data <- tmp %>% 
    abind(along = 4) %>% 
    aperm(c(1,4,2,3))
  
  # --- 2.4. Re-arrange quality checks in an array
  # --- 2.4.1. Colors as values
  qc_col <- mclapply(ensemble_files, function(x){
    memory_cleanup() # low memory use
    val <- x$REC[1:4]
    col <- matrix(data = rep(x$REC$COL, 4), nrow = 6, ncol = 4)
    col[val == 0] <- "white"
    dimnames(col) <- dimnames(val)
    return(col)
  }, mc.cores = MAX_CLUSTERS, mc.cleanup = TRUE) %>% 
    abind(., along = 3) %>% aperm(., c(3,1,2))
  
  # --- 2.4.2. Recommendations
  qc_rec <- mclapply(ensemble_files, function(x){
    memory_cleanup() # low memory use
    rec <- x$REC$Recommandation
    return(rec)
  }, mc.cores = MAX_CLUSTERS, mc.cleanup = TRUE) %>% abind(., along = 2) %>% aperm(., c(2,1))
  dimnames(qc_rec)[[2]] <- dimnames(qc_col)[[2]]
  
  # --- 3. Diversity computing
  # --- 3.1. Identify empty cells
  full_cell_id <- which(!is.na(data[,1,1,1]))
  cell_vector <- terra::values(r0) %>% as.numeric()
  
  # --- 3.2. Compute diversity
  hill_diversity <- mclapply(HILL, function(h){
    memory_cleanup() # low memory use
    tmp <- apply(data, c(3,4), function(x){
      df0 <- x[full_cell_id, ] %>% as.data.frame() # subset full cells
      div0 <- hill_taxa(comm = df0, q = h) # compute diversity
      div <- cell_vector # base
      div[full_cell_id] <- div0 # fill non empty cells
      return(div)
    }) # end apply
  }, mc.cores = MAX_CLUSTERS, mc.cleanup = TRUE) %>% 
    abind(., along = 4) %>% aperm(., c(1,4,2,3))
  
  # add pretty names
  dimnames(hill_diversity)[[2]] <- as.character(HILL)
  dimnames(hill_diversity)[[4]] <- as.character(1:12)
  
  # --- 4. Generate the species NETCDF file
  # --- 4.1. Set up
  library(RNetCDF)
  nc <- create.nc((paste0(project_wd, "/output/", FOLDER_NAME,"/species_output")), format = "netcdf4")
  
  # --- 4.2. Global attributes
  att.put.nc(nc, variable = "NC_GLOBAL", name = "Conventions", type = "NC_CHAR", value = "CF-1.12")
  att.put.nc(nc, variable = "NC_GLOBAL", name = "title", type = "NC_CHAR", value = "Taxa_level_output_CEPHALOPOD")
  att.put.nc(nc, variable = "NC_GLOBAL", name = "institution", type = "NC_CHAR", value = "ETH Zurich")
  att.put.nc(nc, variable = "NC_GLOBAL", name = "source", type = "NC_CHAR", value = "ETH Zurich")
  att.put.nc(nc, variable = "NC_GLOBAL", name = "history", type = "NC_CHAR", value = paste(Sys.time(), "File created"))
  att.put.nc(nc, variable = "NC_GLOBAL", name = "comment", type = "NC_CHAR", value = "Uses attributes recommended by http://cfconventions.org")
  att.put.nc(nc, variable = "NC_GLOBAL", name = "references", type = "NC_CHAR", value = "https://doi.org/10.1111/2041-210X.70040")
  
  # --- 4.3. Dimensions
  # Longitude dimension
  dim.def.nc(nc, dimname = "lon", dimlength = 360)
  var.def.nc(nc, varname = "lon", vartype = "NC_FLOAT", dimensions = "lon")
  att.put.nc(nc, variable = "lon", name = "standard_name", type = "NC_CHAR", value = "longitude")
  att.put.nc(nc, variable = "lon", name = "long_name", type = "NC_CHAR", value = "longitude")
  att.put.nc(nc, variable = "lon", name = "units", type = "NC_CHAR", value = "degrees_east")
  
  # Latitude dimension
  dim.def.nc(nc, dimname = "lat", dimlength = 180)
  var.def.nc(nc, varname = "lat", vartype = "NC_FLOAT", dimensions = "lat")
  att.put.nc(nc, variable = "lat", name = "standard_name", type = "NC_CHAR", value = "latitude")
  att.put.nc(nc, variable = "lat", name = "long_name", type = "NC_CHAR", value = "latitude")
  att.put.nc(nc, variable = "lat", name = "units", type = "NC_CHAR", value = "degrees_north")
  
  # Time dimension
  dim.def.nc(nc, dimname = "time", dimlength = 13)
  var.def.nc(nc, varname = "time", vartype = "NC_FLOAT", dimensions = "time")
  att.put.nc(nc, variable = "time", name = "standard_name", type = "NC_CHAR", value = "month")
  att.put.nc(nc, variable = "time", name = "long_name", type = "NC_CHAR", value = "month")
  att.put.nc(nc, variable = "time", name = "units", type = "NC_CHAR", value = "month from 1 to 12, annual average as 13")
  
  # Target dimension
  dim.def.nc(nc, dimname = "target", dimlength = dim(data)[[2]])
  var.def.nc(nc, varname = "target", vartype = "NC_STRING", dimensions = "target")
  att.put.nc(nc, variable = "target", name = "standard_name", type = "NC_CHAR", value = "target_name")
  att.put.nc(nc, variable = "target", name = "long_name", type = "NC_CHAR", value = "target_name")
  att.put.nc(nc, variable = "target", name = "units", type = "NC_CHAR", value = "target_name")
  
  # Algorithm dimension
  dim.def.nc(nc, dimname = "algorithm", dimlength = dim(qc_col)[[2]])
  var.def.nc(nc, varname = "algorithm", vartype = "NC_STRING", dimensions = "algorithm")
  att.put.nc(nc, variable = "algorithm", name = "standard_name", type = "NC_CHAR", value = "algorithm")
  att.put.nc(nc, variable = "algorithm", name = "long_name", type = "NC_CHAR", value = "algorithm")
  att.put.nc(nc, variable = "algorithm", name = "units", type = "NC_CHAR", value = "algorithm")
  
  # QC dimension
  dim.def.nc(nc, dimname = "qc", dimlength = 4)
  var.def.nc(nc, varname = "qc", vartype = "NC_STRING", dimensions = "qc")
  att.put.nc(nc, variable = "qc", name = "standard_name", type = "NC_CHAR", value = "quality_check")
  att.put.nc(nc, variable = "qc", name = "long_name", type = "NC_CHAR", value = "quality_check")
  att.put.nc(nc, variable = "qc", name = "units", type = "NC_CHAR", value = "quality_check including (1) a priori variable importance, (2) predictive performance, (3) cumulative variance explained, (4) projection uncertainty")
  
  # --- 4.4. Variables
  # CRS variable
  var.def.nc(nc, varname = "crs", vartype = "NC_CHAR", dimensions = NA)
  att.put.nc(nc, variable = "crs", name = "grid_mapping_name", type = "NC_CHAR", value = "latitude_longitude")
  att.put.nc(nc, variable = "crs", name = "long_name", type = "NC_CHAR", value = "CRS definition")
  att.put.nc(nc, variable = "crs", name = "longitude_of_prime_meridian", type = "NC_DOUBLE", value = 0.)
  att.put.nc(nc, variable = "crs", name = "semi_major_axis", type = "NC_DOUBLE", value = 6378137.)
  att.put.nc(nc, variable = "crs", name = "inverse_flattening", type = "NC_DOUBLE", value = 298.257223563)
  att.put.nc(nc, variable = "crs", name = "spatial_ref", type = "NC_CHAR", value = 'GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AXIS[\"Latitude\",NORTH],AXIS[\"Longitude\",EAST],AUTHORITY[\"EPSG\",\"4326\"]]')
  att.put.nc(nc, variable = "crs", name = "GeoTransform", type = "NC_CHAR", value = '-180 0.08333333333333333 0 90 0 -0.08333333333333333 ')
  
  # Annual average projection
  var.def.nc(nc, varname = "mean", vartype = "NC_FLOAT", dimensions = c("lon", "lat", "time", "target"))
  att.put.nc(nc, variable = "mean", name = "long_name", type = "NC_CHAR", value = "Projected annual average by target")
  att.put.nc(nc, variable = "mean", name = "units", type = "NC_CHAR", value = paste("Generated using", CALL$DATA_TYPE, "observations from", CALL$DATA_SOURCE))
  att.put.nc(nc, variable = "mean", name = "grid_mapping", type = "NC_CHAR", value = "crs")
  att.put.nc(nc, variable = "mean", name = "_FillValue", type = "NC_FLOAT", value = -9999.9)
  
  # Standard deviation of the projection
  var.def.nc(nc, varname = "sd", vartype = "NC_FLOAT", dimensions = c("lon", "lat", "time", "target"))
  att.put.nc(nc, variable = "sd", name = "long_name", type = "NC_CHAR", value = "Projected standard deviation across bootstrap replicates by target")
  att.put.nc(nc, variable = "sd", name = "units", type = "NC_CHAR", value = paste("Generated using", CALL$DATA_TYPE, "observations from", CALL$DATA_SOURCE))
  att.put.nc(nc, variable = "sd", name = "grid_mapping", type = "NC_CHAR", value = "crs")
  att.put.nc(nc, variable = "sd", name = "_FillValue", type = "NC_FLOAT", value = -9999.9)
  
  # Quality check colors
  var.def.nc(nc, varname = "qc_col", vartype = "NC_STRING", dimensions = c("target", "algorithm", "qc"))
  att.put.nc(nc, variable = "qc_col", name = "long_name", type = "NC_CHAR", value = "Quality check traffic light by target")
  att.put.nc(nc, variable = "qc_col", name = "units", type = "NC_CHAR", value = "Green, yellow and red show which QC are passed successfully, for 4, 3 or less successfull quality checks")
  
  # Quality check recommendation
  var.def.nc(nc, varname = "qc_rec", vartype = "NC_STRING", dimensions = c("target", "algorithm"))
  att.put.nc(nc, variable = "qc_rec", name = "long_name", type = "NC_CHAR", value = "Quality check recommendationby target")
  att.put.nc(nc, variable = "qc_rec", name = "units", type = "NC_CHAR", value = "Written recommendation based on the traffic light system")
  
  # Target names
  var.def.nc(nc, varname = "target_name", vartype = "NC_STRING", dimensions = c("target"))
  att.put.nc(nc, variable = "target_name", name = "long_name", type = "NC_CHAR", value = "Name of the target")
  
  # Algorithm name
  var.def.nc(nc, varname = "algorithm_name", vartype = "NC_STRING", dimensions = c("algorithm"))
  att.put.nc(nc, variable = "algorithm_name", name = "long_name", type = "NC_CHAR", value = "Name of the algorithm")
  
  # Quality check name
  var.def.nc(nc, varname = "qc_name", vartype = "NC_STRING", dimensions = c("qc"))
  att.put.nc(nc, variable = "qc_name", name = "long_name", type = "NC_CHAR", value = "Name of the quality check")
  
  # Close and check
  sync.nc(nc)
  print.nc(nc)
  
  # --- 4.5. Write values
  # --- 4.5.1. Mean
  proj_mean <- apply(data, c(1,2,4), mean, na.rm = TRUE)
  proj_mean <- abind(proj_mean, apply(data, c(1,2), mean, na.rm = TRUE), along = 3)
  proj_mean <- array(proj_mean, dim = c(360,180,dim(data)[[2]],13)) %>% aperm(., c(1,2,4,3))
  var.put.nc(nc, "mean", proj_mean)
  
  # --- 4.5.2. SD
  proj_sd <- apply(data, c(1,2,4), sd, na.rm = TRUE)
  proj_sd <- abind(proj_sd, apply(proj_sd, c(1,2), mean, na.rm = TRUE), along = 3)
  proj_sd <- array(proj_sd, dim = c(360,180,dim(data)[[2]],13)) %>% aperm(., c(1,2,4,3))
  var.put.nc(nc, "sd", proj_sd)
  
  # --- 4.5.3. Quality checks
  var.put.nc(nc, "qc_col", qc_col)
  var.put.nc(nc, "qc_rec", qc_rec)
  
  # --- 4.5.4. Names
  target_name <- lapply(ensemble_files, function(x){x$SUBFOLDER_NAME}) %>% unlist()
  var.put.nc(nc, "target_name", target_name)
  var.put.nc(nc, "algorithm_name", dimnames(qc_col)[[2]])
  var.put.nc(nc, "qc_name", dimnames(qc_col)[[3]])
  
  # Close and check
  sync.nc(nc)
  print.nc(nc)
  
  close.nc(nc)
  
  # --- 5. Generate the diversity NETCDF file
  # --- 5.1. Set up
  library(RNetCDF)
  nc <- create.nc((paste0(project_wd, "/output/", FOLDER_NAME,"/diversity_output")), format = "netcdf4")
  
  # --- 5.2. Global attributes
  att.put.nc(nc, variable = "NC_GLOBAL", name = "Conventions", type = "NC_CHAR", value = "CF-1.12")
  att.put.nc(nc, variable = "NC_GLOBAL", name = "title", type = "NC_CHAR", value = "Divrsity_output_CEPHALOPOD")
  att.put.nc(nc, variable = "NC_GLOBAL", name = "institution", type = "NC_CHAR", value = "ETH Zurich")
  att.put.nc(nc, variable = "NC_GLOBAL", name = "source", type = "NC_CHAR", value = "ETH Zurich")
  att.put.nc(nc, variable = "NC_GLOBAL", name = "history", type = "NC_CHAR", value = paste(Sys.time(), "File created"))
  att.put.nc(nc, variable = "NC_GLOBAL", name = "comment", type = "NC_CHAR", value = "Uses attributes recommended by http://cfconventions.org")
  att.put.nc(nc, variable = "NC_GLOBAL", name = "references", type = "NC_CHAR", value = "https://doi.org/10.1111/2041-210X.70040")
  
  # --- 5.3. Dimensions
  # Longitude dimension
  dim.def.nc(nc, dimname = "lon", dimlength = 360)
  var.def.nc(nc, varname = "lon", vartype = "NC_FLOAT", dimensions = "lon")
  att.put.nc(nc, variable = "lon", name = "standard_name", type = "NC_CHAR", value = "longitude")
  att.put.nc(nc, variable = "lon", name = "long_name", type = "NC_CHAR", value = "longitude")
  att.put.nc(nc, variable = "lon", name = "units", type = "NC_CHAR", value = "degrees_east")
  
  # Latitude dimension
  dim.def.nc(nc, dimname = "lat", dimlength = 180)
  var.def.nc(nc, varname = "lat", vartype = "NC_FLOAT", dimensions = "lat")
  att.put.nc(nc, variable = "lat", name = "standard_name", type = "NC_CHAR", value = "latitude")
  att.put.nc(nc, variable = "lat", name = "long_name", type = "NC_CHAR", value = "latitude")
  att.put.nc(nc, variable = "lat", name = "units", type = "NC_CHAR", value = "degrees_north")
  
  # Time dimension
  dim.def.nc(nc, dimname = "time", dimlength = 13)
  var.def.nc(nc, varname = "time", vartype = "NC_FLOAT", dimensions = "time")
  att.put.nc(nc, variable = "time", name = "standard_name", type = "NC_CHAR", value = "month")
  att.put.nc(nc, variable = "time", name = "long_name", type = "NC_CHAR", value = "month")
  att.put.nc(nc, variable = "time", name = "units", type = "NC_CHAR", value = "month from 1 to 12, annual average as 13")
  
  # Target dimension
  dim.def.nc(nc, dimname = "target", dimlength = dim(hill_diversity)[[2]])
  var.def.nc(nc, varname = "target", vartype = "NC_STRING", dimensions = "target")
  att.put.nc(nc, variable = "target", name = "standard_name", type = "NC_CHAR", value = "target_name")
  att.put.nc(nc, variable = "target", name = "long_name", type = "NC_CHAR", value = "target_name")
  att.put.nc(nc, variable = "target", name = "units", type = "NC_CHAR", value = "target_name")
  
  # --- 5.4. Variables
  # CRS variable
  var.def.nc(nc, varname = "crs", vartype = "NC_CHAR", dimensions = NA)
  att.put.nc(nc, variable = "crs", name = "grid_mapping_name", type = "NC_CHAR", value = "latitude_longitude")
  att.put.nc(nc, variable = "crs", name = "long_name", type = "NC_CHAR", value = "CRS definition")
  att.put.nc(nc, variable = "crs", name = "longitude_of_prime_meridian", type = "NC_DOUBLE", value = 0.)
  att.put.nc(nc, variable = "crs", name = "semi_major_axis", type = "NC_DOUBLE", value = 6378137.)
  att.put.nc(nc, variable = "crs", name = "inverse_flattening", type = "NC_DOUBLE", value = 298.257223563)
  att.put.nc(nc, variable = "crs", name = "spatial_ref", type = "NC_CHAR", value = 'GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AXIS[\"Latitude\",NORTH],AXIS[\"Longitude\",EAST],AUTHORITY[\"EPSG\",\"4326\"]]')
  att.put.nc(nc, variable = "crs", name = "GeoTransform", type = "NC_CHAR", value = '-180 0.08333333333333333 0 90 0 -0.08333333333333333 ')
  
  # Annual average projection
  var.def.nc(nc, varname = "mean", vartype = "NC_FLOAT", dimensions = c("lon", "lat", "time", "target"))
  att.put.nc(nc, variable = "mean", name = "long_name", type = "NC_CHAR", value = "Projected annual average by target")
  att.put.nc(nc, variable = "mean", name = "units", type = "NC_CHAR", value = paste("Generated using", CALL$DATA_TYPE, "observations from", CALL$DATA_SOURCE))
  att.put.nc(nc, variable = "mean", name = "grid_mapping", type = "NC_CHAR", value = "crs")
  att.put.nc(nc, variable = "mean", name = "_FillValue", type = "NC_FLOAT", value = -9999.9)
  
  # Standard deviation of the projection
  var.def.nc(nc, varname = "sd", vartype = "NC_FLOAT", dimensions = c("lon", "lat", "time", "target"))
  att.put.nc(nc, variable = "sd", name = "long_name", type = "NC_CHAR", value = "Projected standard deviation across bootstrap replicates by target")
  att.put.nc(nc, variable = "sd", name = "units", type = "NC_CHAR", value = paste("Generated using", CALL$DATA_TYPE, "observations from", CALL$DATA_SOURCE))
  att.put.nc(nc, variable = "sd", name = "grid_mapping", type = "NC_CHAR", value = "crs")
  att.put.nc(nc, variable = "sd", name = "_FillValue", type = "NC_FLOAT", value = -9999.9)
  
  # Target names
  var.def.nc(nc, varname = "target_name", vartype = "NC_STRING", dimensions = c("target"))
  att.put.nc(nc, variable = "target_name", name = "long_name", type = "NC_CHAR", value = "Name of the target")
  
  # Close and check
  sync.nc(nc)
  print.nc(nc)
  
  # --- 5.5. Write values
  # --- 5.5.1. Mean
  proj_mean <- apply(hill_diversity, c(1,2,4), mean, na.rm = TRUE)
  proj_mean <- abind(proj_mean, apply(hill_diversity, c(1,2), mean, na.rm = TRUE), along = 3)
  proj_mean <- array(proj_mean, dim = c(360,180,dim(hill_diversity)[[2]],13)) %>% aperm(., c(1,2,4,3))
  var.put.nc(nc, "mean", proj_mean)
  
  # --- 5.5.2. SD
  proj_sd <- apply(hill_diversity, c(1,2,4), sd, na.rm = TRUE)
  proj_sd <- abind(proj_sd, apply(proj_sd, c(1,2), mean, na.rm = TRUE), along = 3)
  proj_sd <- array(proj_sd, dim = c(360,180,dim(hill_diversity)[[2]],13)) %>% aperm(., c(1,2,4,3))
  var.put.nc(nc, "sd", proj_sd)
  
  # --- 5.5.4. Names
  target_name <- paste0("Hill_scaling_factor_", HILL)
  var.put.nc(nc, "target_name", target_name)

  # Close and check
  sync.nc(nc)
  print.nc(nc)
  
  close.nc(nc)
  
  # --- 6. Plot diversities
  pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/06_hill_diversity.pdf"))

  for(h in seq_along(HILL)){
    for(m in seq_along(MONTH)){
      # --- 3.1. Format the data
      div <- hill_diversity[,HILL[h],,MONTH[[m]]]
      div <- apply(div, 1, mean, na.rm = TRUE)
      
      # --- 3.3. Plot
      r <- terra::setValues(r0, div)
      plot(r, col = parula_pal(100), main = paste("Hill diversity of", HILL[h], "(nb. of species)", "\n Month:", paste(MONTH[[m]], collapse = ",")))
      plot(land, col = "black", add = TRUE, legend = FALSE)
      
    } # end month loop
  } # end hill loop
  
  dev.off()
  
} # END FUNCTION


