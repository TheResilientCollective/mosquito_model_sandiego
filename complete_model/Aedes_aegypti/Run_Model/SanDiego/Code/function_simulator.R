## Function simulator
## Arguments:
## 1. initial conditions
## 2. LP_mortality alpha1
## 3. LP_mortality beta1
## 4. Oviposition

mosq_sim <- function(alpha1_val = -2.514493709671852,
                     beta1_val = -4.457398864610944,
                     init_val = 0.01,
                     lifestage = "ActMosq", # Stage options "ActMosq", "TotMosq", "LP", "Eggs"
                     config_file = "/Users/hwaddel/Documents/Research/Mosquito_Modeling/Code/mosquito_model_sandiego/complete_model/Aedes_aegypti/Run_Model/SanDiego/simulator_config.yaml"){ 
  
  options(warn = -1)
  suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
  
  library(yaml,quietly = TRUE)
  options("rgdal_show_exportToProj4_warnings"="none")
  library(sp,quietly = TRUE)
  library(logger,quietly = TRUE) 
  library(tidyverse,quietly = TRUE)
  library(lubridate,quietly = TRUE)
  library(stringr,quietly = TRUE)
  library(ggplot2,quietly = TRUE)
  library(ggpubr,quietly = TRUE)
  library(dplyr,quietly = TRUE)
  
  ### Read in Config File ###
  config_data <- yaml.load_file(config_file, eval.expr = TRUE)

  source(paste0(config_data$MODEL_BASE_PATH,"Code/mosquito_model_source/all_run_fxns.R"), chdir = T)
  source(paste0(config_data$MODEL_BASE_PATH,"Code/mosquito_model_source/mosq_param_funcs.R"), chdir = T)
  source(paste0(config_data$MODEL_BASE_PATH,"Code/mosquito_model_source/data_import_function.R"), chdir = T)
  
  mosq_type = config_data$MOSQ_TYPE
  pop_type = config_data$POP_TYPE
  burnin = config_data$BURNIN
  
  if (mosq_type == 'culex'){
    # parameters that are constant across runs
    params <- config_data$PARAMETERS$CULEX$MOSQ_PARAMS
    
    # initial values for parameters that we fit 
    fit_params <- config_data$PARAMETERS$CULEX$FIT_PARAMS
    
    # numerical solver
    source(paste0(config_data$MODEL_BASE_PATH,"/Code/mosquito_model_source/culex_ss_test.R"))
  }
  
  if (mosq_type == 'a_aegypti'){
    # parameters that are constant across runs
    params <- config_data$PARAMETERS$A_AEGYPTI$MOSQ_PARAMS
    
    # initial values for parameters that we fit
    fit_params <- config_data$PARAMETERS$A_AEGYPTI$FIT_PARAMS
    
    # numerical solver
    source(paste0(config_data$MODEL_BASE_PATH,"/Code/mosquito_model_source/aegypti_ss_test.R"))
  }
  
  N <- config_data$N_VAL  # number of age bins
  
  config_data$MOSQ_DATA_PATH <- paste0(config_data$MODEL_BASE_PATH, config_data$MOSQ_DATA_PATH)
  config_data$EXT_DATA_PATH <- paste0(config_data$MODEL_BASE_PATH, config_data$EXT_DATA_PATH)
  
  input.data <- data.import(config_data) # read in config file
  observed.data <- input.data$observed  # trap data
  external.signal <- input.data$external  # temp, water data, and daylight hours
  
  # evap_etc <- read.csv("~/Documents/Research/Mosquito_Modeling/Code/culex-and-aedes-aegypti/input-data/SanDiego/SanDiego_PBM_input_evaporation_etc_0_2015to2026.csv")
  # other_wetness_var <- evap_etc$EvapMean[evap_etc$date %in% external.signal$Date]
  
  # external.signal$Wetness <- external.signal$PrecipMean / (max(external.signal$PrecipMean) - min(external.signal$PrecipMean))
  # external.signal$Wetness <- (other_wetness_var - min(other_wetness_var))/(max(other_wetness_var) - min(other_wetness_var))
  
  if(as.numeric(config_data$TEMP_LAG) > 0){
    temp_lag <- as.numeric(config_data$TEMP_LAG)
    ## Adjust the column up by x
    Temp_tmp <- external.signal$Temperature[1:(length(external.signal$Temperature) - temp_lag)]
    Temp_tmp <- c(external.signal$Temperature[1:temp_lag],Temp_tmp)
    external.signal$Temperature <- Temp_tmp
  }
  
  if(as.numeric(config_data$WETNESS$lag) > 0){
    wet_lag <- as.numeric(config_data$WETNESS$lag)
    ## Adjust the column up by x
    Wetness_tmp <- external.signal$Wetness[1:(length(external.signal$Wetness) - wet_lag)]
    Wetness_tmp <- c(external.signal$Wetness[1:wet_lag],Wetness_tmp)
    external.signal$Wetness <- Wetness_tmp
  }
  
  config_data$MODEL_RUN_TYPE <- 'test'
  
  # read in fitted params file
  fitted_params <- read_csv(paste0(config_data$MODEL_BASE_PATH,config_data$FIT_PARAM_DATA_PATH), show_col_types = F)
  
  # we want the fitted_params content as a named vector
  fitted_param_vals <- fitted_params$fit_parameters
  names(fitted_param_vals) <- fitted_params$parameter_names
  
  fitted_param_vals["LP_Mortality.alpha1"] <- alpha1_val
  fitted_param_vals["LP_Mortality.beta1"] <- beta1_val
  
  # function found in all_run_fxns.R
  out <- run_model_input(fitted_param_vals, params, external.signal, N)
  
  mosq_vec <- out %>%
    filter(times >= burnin*365) %>%
    mutate(lMosq = log(get(lifestage) + 0.5)) %>%
    dplyr::select(lMosq) %>%
    unlist(use.names = F)
  
  return(mosq_vec)
  # cat(paste(actmosq_vec, collapse = " "), "\n")
  # cat(jsonlite::toJSON(actmosq_vec,dights=NA))
}
