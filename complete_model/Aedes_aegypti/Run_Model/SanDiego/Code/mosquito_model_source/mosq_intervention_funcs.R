library(tidyverse)

calc_exponential_intervention <- function(daystep, application_days, max_percent, decay_exponent, max_duration){
  if(max_percent < 0 | max_percent > 1){
    
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
  intervention_progress <- daystep - app_day
  
  if(intervention_progress < 0){
    stop("Daystep error, intervention")
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
    return(max_percent)
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
    return(calc_exponential_intervention(daystep_arg, application_days_arg, percent_effective_arg, decay_exponent_arg, max_duration_arg))
  } else if(intervention_function_type == "linear"){
    return(calc_linear_intervention(daystep_arg, application_days_arg, percent_effective_arg, max_duration_arg))
  } else{
    stop("unknown intervention function")
  }
}
