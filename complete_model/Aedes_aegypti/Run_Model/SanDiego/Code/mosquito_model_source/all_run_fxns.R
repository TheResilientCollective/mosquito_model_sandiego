## Source Model Functions
library(subplex)
# library(optimr)
library(optimx)

source("./culex_simulation_function.R")
source("./convolution_function.R")

# this is used in run_model_input and fig_obj functions below
run_model <- function(params, external.signal, N){
  input1_Temp <- external.signal$Temperature
  input2_DayHrs <- external.signal$DaylightHrs
  input3_Wet <- external.signal$Wetness
  
  times <- external.signal$TimeSeq
  daysteps <- external.signal$Daystep
  
  model_df <- CulexSimFunc(times, params, input1_Temp, 
                           input2_DayHrs, input3_Wet, N)
  
  model_out <- cbind(data.frame(model_df), times, 
                     input1_Temp, input2_DayHrs, 
                     input3_Wet, daysteps)
  return(model_out)
}

# this runs the run_model function defined above
# found in models_main.R for testing
run_model_input <- function(input_fitted_params, params, external.signal, N){
  
  # combine fitted params with constant params
  unlist_params <- unlist(params)
  unlist_params[names(input_fitted_params)] <- input_fitted_params
  new_params <- relist(unlist_params, params)
  
  model_out <- run_model(new_params, external.signal, N)
  return(model_out)
}

# this also runs the run_model function defined above
# found in fit_model function below
fit_obj <- function(input_fitted_params, params, external.signal, 
                    observed.data, N, pop_type){
  
  # add initial guesses of fitted parameters into the parameter list
  unlist_params <- unlist(params)
  unlist_params[names(input_fitted_params)] <- input_fitted_params
  new_params <- relist(unlist_params, params)
  
  out <- run_model(new_params, external.signal, N)
  
  mosq_out <- out$ActMosq  # model output
  
  daysteps <- observed.data$Daystep
  ave_mosq_count <- observed.data$muMosq  # trap data
  cost_weight <- observed.data$cost_weight
  
  # shift output so that 1st observed data point aligns with model output
  # will work as long as external data starts before or on the 1st observed data point
  shift <- which(external.signal$Daystep == min(daysteps)) 
  Y = mosq_out[c(daysteps - min(daysteps) + shift)]
  
  Y[is.na(Y) == TRUE] = 0.0    # set NA values to 0
  
  if(pop_type == 'raw'){
    
    Y[Y <= 0] = 0 # threshold for non-logged populations
    SSq = sum((cost_weight*(Y - ave_mosq_count))^2) # ave_mosq_count needs to be the non-logged data
    
  } else if (pop_type == 'log'){
    
    Y[Y <= 0.0] = 0.5   # threshold for logged populations
    SSq = sum(cost_weight*(log(Y) - ave_mosq_count)^2) # ave_mosq_count needs to be the logged data
  }
  return(SSq)
}

# found in models_main.R for training
# runs parameter optimization
fit_model <- function(input_fit_params, params, external.signal, observed.data, N, pop_type){
  
  check <- fit_obj(unlist(input_fit_params), params, external.signal, observed.data, N, pop_type)

  # optimiztation algorithm
  if(is.numeric(check)){
    print("Passed Check, starting full fit process.")
    
    ptm = proc.time()
    set.seed(45) 
    
    Fit1 = optim(par = unlist(input_fit_params),
                 fn = fit_obj,
                 params = params,
                 external.signal = external.signal,
                 observed.data = observed.data,
                 N = N,
                 pop_type = pop_type,
                 method = "Nelder-Mead",
                 control = list(maxit = 10000), hessian = TRUE)
    
    time_elapsed <- as.vector(proc.time() - ptm)
    cat("Fit Time: ", round(time_elapsed[3], 4), "\n")
    
    if(Fit1$convergence == 0){
      
      param_fits <- c(Fit1$par)
      return(param_fits)
      
    } else {
      print("Parameter Fit Function DID NOT CONVERGE")
    }
    
  } else {
    print("SSqFunc did not produce a numerical result")
  }
}
