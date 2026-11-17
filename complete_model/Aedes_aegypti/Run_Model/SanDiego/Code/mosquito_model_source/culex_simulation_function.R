## Culex Abundance Simulation 

### Arguments ###
# 1: times - timesteps
# 2: parms - list of model parameters
# 3: import1_Temp - vector of mean daily temperature
# 4: import2_DayHrs - vector of normalized daily daylight hrs
# 5: import3_Wet - vector of mean normalized water data

### Final_Output ###
# returns populations for Active mosquitoes, Total mosquitoes, eggs, and juveniles

CulexSimFunc = function(times, parms, import1_Temp, import2_DayHrs, import3_Wet, N = 100){
  init = as.numeric(parms['init'][[1]])
  
  force_init_positive = F

  if(force_init_positive){
    init = exp(as.numeric(parms['init'][[1]]))
  }
  
  # mosq_type = "a_aegypti"
  
  # Initialize results matrix
  if(length(parms$Diapause) == 0){
    mosq_type = "a_aegypti"
    # OutMat = matrix(rep(0.0, 2*length(times)), nrow = length(times), ncol = 19)
    OutMat = matrix(data=0.0, nrow = length(times), ncol = 21)
    
    colnames(OutMat) = c('ActMosq', 'TotMosq','Eggs','LP','AD1','AD2','AD3','AD4','AddEggs','AddLP','AddAD1','AddAD2','AddAD3','AddAD4','MortLP','MortAD','Ovipos','NOvipos','vtE','vtLP','vtA')
  } else{
    mosq_type = "culex"
    OutMat = matrix(rep(0.0, 2*length(times)), nrow = length(times), ncol = 4)
    
    colnames(OutMat) = c('ActMosq', 'TotMosq','Eggs','LP')
  }
  
  # Initial conditions:
  yinout = rep(0.0, 6*N)
  yinout[(2*N + 1):(6*N)] = init
  
  for(i in 1:length(times)){
    out = cpmod(time = 1.0, daystep = times[i], yinout = yinout , parms = parms, 
                   import1_Temp = import1_Temp[i], import2_DayHrs = import2_DayHrs[i], 
                   import3_Wet = import3_Wet[i], N = N)
    
    yinout = out$yinout
    
    # Updating the results matrix
    OutMat[i,'ActMosq']  <- out$data[['ActMosq']][[1]]
    OutMat[i,'TotMosq']  <- out$data[['TotMosq']][[1]]
    OutMat[i,'Eggs']  <- out$data[['Eggs']][[1]]
    OutMat[i,'LP']  <- out$data[['LP']][[1]]
    
    if(mosq_type == "a_aegypti"){
      OutMat[i,'AD1'] <- out$data[['AD1']][[1]]
      OutMat[i,'AD2'] <- out$data[['AD2']][[1]]
      OutMat[i,'AD3'] <- out$data[['AD3']][[1]]
      OutMat[i,'AD4'] <- out$data[['AD4']][[1]]
      OutMat[i,'AddEggs']  <- out$data[['AddEggs']][[1]]
      OutMat[i,'AddLP']  <- out$data[['AddLP']][[1]]
      OutMat[i,'AddAD1']  <- out$data[['AddAD1']][[1]]
      OutMat[i,'AddAD2']  <- out$data[['AddAD2']][[1]]
      OutMat[i,'AddAD3']  <- out$data[['AddAD3']][[1]]
      OutMat[i,'AddAD4']  <- out$data[['AddAD4']][[1]]
      OutMat[i,'MortLP']  <- out$data[['MortLP']][[1]]
      OutMat[i,'MortAD']  <- out$data[['MortAD']][[1]]
      OutMat[i,'Ovipos']  <- out$data[['Ovipos']][[1]]
      OutMat[i,'NOvipos']  <- out$data[['NOvipos']][[1]]
      OutMat[i,'vtE'] <- out$data[['vtE']][[1]]
      OutMat[i,'vtLP'] <- out$data[['vtLP']][[1]]
      OutMat[i,'vtA'] <- out$data[['vtA']][[1]]
    }
  }

  return(OutMat)
}
