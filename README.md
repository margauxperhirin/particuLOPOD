## **particuLOPOD - Adaped pipeline of CEPHALOPOD for multi-outputs**

This branch adds the possibly to predict a multi-output (vector of 3 values) in the CEPHALOPOD habitat modelling pipeline, developed for the Ecosystem - Workbench of the Bluecloud2026 E.U. project by Schickele et al. 2025.**

This wad developped for Perhirin et al. (in prep.) to predict 3 coefficients together in the same models. 

**Main changes**

* possibility to predict 3 outputs from one model

* normalisation of environmental data before running the MPL

**Important points**

* only the continuous pipeline was modified

* RF and MLP are the ONLY models that can be run in particuLOPOD

* your dataset must now contain the following variables:
  
        - worms_id : AphiaID or other identifier for the species/taxon
  
        - decimallatitude : latitude of the sample in decimal degrees (-90 to +90)
  
        - decimallongitude : longitude of the sample in decimal degrees (-180 to +180)
  
        - depth : sample depth in meters
  
        - year : year of sampling (integer)
  
        - month : month of sampling (integer)
  
        - measurementunit : units of the measurement value
  
        - taxonrank : taxonomic rank (e.g., species, genus, order...)
  
        - PSD1 : first coefficient to predict
  
        - PSD2 : second coefficient to predict
  
        - PSD3 : third coefficient to predict  


**What was modified?**

* 00_config: addition of the package *randomforestSRC* 
* 01d_list_custom: psd1, psd2 and psd3 in the required columns
* 03d_query_custom: psd1, psd2 and psd3 in Y but removed in S
* 04_query_env: creation of normalised predictors and their parameters for MLP (stacked in QUERY under X_norm and X_norm_params)
* 05_pseudo_abs: psd1, psd2 and psd3 instead of one variable for the plots in 01_observations.pdf but no modification concerning pseudo-absences
* 09b_model_continuous: addition of specific cases when model = RF or model = MLP
* 10b_eval_continuous: addition of specific cases when model = RF or model = MLP, computation of R2 and RMSE standard deviations and means (among the 3 coefficients)
* 11b_proj_continuous: addition of specific cases when model = RF or model = MLP for the loops on bootstraps and on months
* 12a_standardmaps: supplementary loop on the 3 coefficients


Last modification on the 6th of May 2026.
