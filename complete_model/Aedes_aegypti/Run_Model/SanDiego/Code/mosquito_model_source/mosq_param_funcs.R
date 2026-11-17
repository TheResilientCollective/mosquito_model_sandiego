library(tidyverse)

eyring <- function(temp, params){
  return(as.numeric(as.numeric(params['Psi'][[1]])*(273.15 + temp)*exp(-as.numeric(params['AE'][[1]])/(8.3144598*(273.15+temp)))))
}

briere <- function(temp, params){
  return(as.numeric(as.numeric(params['c'][[1]])*temp*(temp - as.numeric(params['T0'][[1]]))*sqrt(as.numeric(params['Tm'][[1]]) - temp)))
}

briere_thresholded <- function(temp, params){
  if(temp >= as.numeric(params['Tm'][[1]])){
    return(0.01)
  } else{
    return(as.numeric(as.numeric(params['c'][[1]])*temp*(temp - as.numeric(params['T0'][[1]]))*sqrt(as.numeric(params['Tm'][[1]]) - temp))) 
  }
}

linear <- function(temp, params){
  return(as.numeric(as.numeric(params['Int'][[1]]) + as.numeric(params['Coef'][[1]])*temp))
}

sharpe_schoolfield <- function(temp, params){
  return((as.numeric(params['r'][[1]]) * ((temp + 273.15)/298) * exp((as.numeric(params['Ha'][[1]])/1.987) * (1/298 - 1/(temp + 273.15)))) / (1 + exp((as.numeric(params['Hh'][[1]])/1.987) * (1/as.numeric(params['T1_2'][[1]]) - 1/(temp + 273.15)))))
}

# 6-parameter version of sharpe-schoolfield
sharpe_schoolfield_6 <- function(temp, params){
  return( ( as.numeric(params['p_25c'][[1]]) * (temp + 273.15)/298 * exp( ( as.numeric(params['H_A'][[1]])/1.987) * (1/298 - 1/(temp + 273.15)) ))/(1 + exp( as.numeric(params['H_L'][[1]])/1.987 * (1/ as.numeric(params['T_12L'][[1]]) - 1/(temp+273.15))) + exp(as.numeric(params['H_H'][[1]])/1.987 * (1/as.numeric(params['T_12H'][[1]]) - 1/(temp + 273.15)))) + 0.025 )
}

exponential <- function(wet, params){
  return(exp(as.numeric(params['alpha1'][[1]]) + as.numeric(params['beta1'][[1]])*wet))
}

diapause_exp <- function(dayhrs, params){
  return(1.0/(1.0 + exp(-(as.numeric(params['Int'][[1]]) - as.numeric(params['Coef'][[1]])*dayhrs))))
}

## helper to find quadratic roots

quadRoots <- function(params) {
  a <- as.numeric(params['a'][[1]])
  b <- as.numeric(params['b'][[1]])
  c <- as.numeric(params['c'][[1]])
  
  discriminant <- (b^2) - (4*a*c)
  
  if(discriminant < 0) {
    return(NULL)  # No real roots
  } else if(discriminant > 0) {
    x1 <- (-b + sqrt(discriminant)) / (2*a)
    x2 <- (-b - sqrt(discriminant)) / (2*a)
    return(c(x1, x2))
  } else {
    x <- -b / (2*a)
    return(c(x, x))  # Both roots are the same
  }
}

constant <- function(temp, params){
  return(as.numeric(params['val'][[1]]))
}

exp_quadratic <- function(temp, peak = 0.41, t_opt = 27, 
                          sd_left = 3.3, sd_right = 3.8){
  return(params['peak'] * exp(-0.5 * ((temp - params['t_opt']) / params['sd'])^2))
}

# quadratic: y = ax^2 + bx + c
quadratic <- function(temp, params){
  roots <- quadRoots(params)
  # print(quadRoots)
  
  # If roots is NULL (no real roots), return 0
  if (is.null(roots)) {
    return(0)
  }
  
  if (temp < roots[1] || temp > roots[2]) {
    return(0)
  } else {
    a <- as.numeric(params['a'][[1]])
    b <- as.numeric(params['b'][[1]])
    c <- as.numeric(params['c'][[1]])
    return(a * temp^2 + b * temp + c)
  }
}

quartic <- function(temp, params){
  b0 = as.numeric(params['b0'][[1]])
  b1 = as.numeric(params['b1'][[1]])
  b2 = as.numeric(params['b2'][[1]])
  b3 = as.numeric(params['b3'][[1]])
  b4 = as.numeric(params['b4'][[1]])
  
  if (temp < 11.623){
    return(0)
  } else {
    return (b0 + b1 * temp + b2 * temp^2 + b3 * temp^3 + b4 * temp^4) }
}

### Sigmoidal function for oviposition
sigmoidal <- function(x, min.val, max.val){
  logistic = 1/(1 + exp(-x))
  
  return(min.val + (logistic*(max.val - min.val)))
}

calc_mosq_atr <- function(env_var, params){
  if (params['func'][[1]] == 'eyring'){
    val <- eyring(env_var, params)
  }
  if (params['func'][[1]] == 'briere'){
    val <- briere(env_var, params)
  }
  if (params['func'][[1]] == 'briere_thresholded'){
    val <- briere_thresholded(env_var, params)
  }
  if (params['func'][[1]] == 'linear'){
    val <- linear(env_var, params)
  }
  if (params['func'][[1]] == 'sharpe-schoolfield'){
    val <- sharpe_schoolfield(env_var, params)
  }
  if (params['func'][[1]] == 'sharpe-schoolfield-6'){
    val <- sharpe_schoolfield_6(env_var, params)
  }
  if (params['func'][[1]] == 'exponential'){
    val <- exponential(env_var, params)
  }
  if (params['func'][[1]] == 'diapause_exp'){
    val <- diapause_exp(env_var, params)
  }
  if (params['func'][[1]] == 'quadratic') {
    val <- quadratic(env_var, params)
  }
  if (params['func'][[1]] == 'quartic') {
    val <- quartic(env_var, params)
  }
  if (params['func'][[1]] == 'constant'){
    val <- constant(env_var, params)
  }
  return(val)
}

###### INTERVENTION FUNCTIONS:
### FUNCTIONS THAT CALCULATE INTERVENTION EFFECTIVENESS
calc_exponential_intervention <- function(daystep, application_days, max_percent, decay_exponent, max_duration){
  if(max_percent < 0 | max_percent > 1){
    stop('Max Percent must be between 0 and 1')
  }
  
  if(decay_exponent < 0){
    stop("Decay exponent must be > 0")
  }
  
  # Intervention not applied yet
  if(daystep < min(application_days)){
    return(0.0)
  }
  
  # find closest application day prior to current daystep
  app_day <- max(application_days[(application_days - daystep <= 0)])
  
  # If duration has ended, no intervention
  if(app_day + max_duration < daystep){
    return(0.0)
  }
  
  # Number of days intervention has been in effect
  intervention_progress <- daystep - app_day + 1
  
  if(intervention_progress <= 0){
    return(0.0)
  }
  
  return(max_percent * (intervention_progress^(-decay_exponent)))
}

calc_linear_intervention <- function(daystep, application_days, max_percent, max_duration){
  # Intervention not applied yet
  if(daystep < min(application_days)){
    return(0.0)
  }
  
  # find closest application day prior to current daystep
  app_day <- max(application_days[(application_days - daystep <= 0)])
  
  # If duration has ended, no intervention
  if(app_day + max_duration < daystep){
    return(0.0)
  }
  
  # Number of days intervention has been in effect
  intervention_progress <- daystep - app_day
  
  if(intervention_progress < 0){
    stop("Daystep error, intervention")
    return(0.0)
  }
  
  if(intervention_progress > max_duration){
    return(0.0)
  }
  
  return(max_percent * (max_duration - intervention_progress)/max_duration)
}

calc_constant_intervention <- function(daystep, application_days, percent_effective){
  if(daystep %in% application_days){
    return(percent_effective)
  } else{
    return(0.0)
  }
}

calc_custom_intervention <- function(daystep, application_days, percent_effective){
  if(!(daystep %in% application_days)){
    return(0.0)
  }
  if(length(application_days) != length(percent_effective)){
    stop("days and effectiveness must have same length")
  }
  
  if(length(unique(application_days)) < length(application_days)){
    stop("provide list of unique days")
  }
  
  if(max(percent_effective) < 0 | max(percent_effective) > 1){
    stop("percent effective must be between 0 and 1")
  }
  
  return(percent_effective[which(application_days == daystep)])
}

calc_intervention_percent <- function(daystep_arg, intervention_function_type, application_days_arg, percent_effective_arg, max_duration_arg, decay_exponent_arg){
  if(intervention_function_type == "constant"){
    return(calc_constant_intervention(daystep_arg, application_days_arg, percent_effective_arg))
  } else if(intervention_function_type == "custom"){
    return(calc_custom_intervention(daystep_arg, application_days_arg, percent_effective_arg))
  } else if(intervention_function_type == "exponential"){
    return(calc_exponential_intervention(daystep = daystep_arg, application_days = application_days_arg, max_percent = percent_effective_arg, decay_exponent = decay_exponent_arg, max_duration = max_duration_arg))
  } else if(intervention_function_type == "linear"){
    return(calc_linear_intervention(daystep_arg, application_days_arg, percent_effective_arg, max_duration_arg))
  } else{
    stop("unknown intervention function")
  }
}
