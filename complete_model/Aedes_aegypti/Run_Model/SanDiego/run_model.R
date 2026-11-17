## Script to run one iteration of the mosquito model
library(yaml)

args <- commandArgs(trailingOnly = T)

if(length(args) == 0){
  alpha1_val = -2.514493709671852
  beta1_val = -4.457398864610944
  init_val = 0.01
  lifestage = "ActMosq" # Stage options "ActMosq", "TotMosq", "LP", "Eggs"
  config_file = "/Users/hwaddel/Documents/Research/Mosquito_Modeling/Code/mosquito_model_sandiego/complete_model/Aedes_aegypti/Run_Model/SanDiego/simulator_config.yaml"
} else if(length(args) == 5){
  alpha1_val = as.numeric(args[1])
  beta1_val = as.numeric(args[2])
  init_val = as.numeric(args[3])
  lifestage = as.character(args[4])
  config_file = as.character(args[5])

  if(!(lifestage %in% c("ActMosq", "TotMosq", "LP", "Eggs"))){
    stop("Unknown lifestage")
  }
} else{
  stop("Wrong number of arguments")
}

config_data <- yaml.load_file(config_file, eval.expr = TRUE)

setwd(config_data$MODEL_BASE_PATH)

source("./Code/function_simulator.R", chdir = T)

est_mosq_pop <- mosq_sim(alpha1_val, beta1_val, init_val, lifestage, config_file)

cat(paste0(est_mosq_pop, collapse = "\t"))