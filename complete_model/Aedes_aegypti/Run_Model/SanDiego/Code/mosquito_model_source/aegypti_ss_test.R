## Numerical method for culex populations

### Arguments: ###
# 1: time - timesteps
# 2: yinout - vector of populations for each discrete age (development) bin
# concatenates discrete age (development) bins for egg, juvenile, 
# and the four adult life stages
# 3: parms - list of model parameters
# 4: import1_Temp - vector of mean daily temperature
# 5: import2_DayHrs - vector of normalized daily daylight hrs
# 6: import3_Wet - vector of mean normalized water data
# 7: N - total number of development bins per life stage

### Final Output: ###
#  yinout: population of each development bin for the next run
#  data: active and total mosquito populations

cpmod = function(time, daystep, yinout, parms, import1_Temp, import2_DayHrs, import3_Wet, N = 100){
  # Setting up the domain: number of bins in the discretized domain for each life stage
  da = 1.0/N  # bin thickness
  
  # kill off juvenile stages and eggs if temp < 0
  # if(import1_Temp < 0.0){
  #   yinout[1:(2*N)] = 0.0
  # }
  
  if(import1_Temp < -15.0){
    yinout[1:(2*N)] = 0.0
  }
  
  # subsetting the different life stages
  Eggs = yinout[1:(N)] # Egg stage
  LP = yinout[(N + 1):(2*N)] # L1 to pupa stages
  AD1 = yinout[(2*N + 1):(3*N)] # first gonotrophic cycle of adult stage
  AD2 = yinout[(3*N + 1):(4*N)] # second gonotrophic cycle of adult stage
  AD3 = yinout[(4*N + 1):(5*N)] # third gonotrophic cycle of adult stage
  AD4 = yinout[(5*N + 1):(6*N)] # fourth gonotrophic cycle of adult stage
  
  # development bins for each life stage
  ageEggs = seq(da, 1.0, by = da)
  ageLP = seq(da, 1.0, by = da)
  ageAD1 = seq(da, 1.0, by = da)
  
  AddEggs = 0.0
  AddLP = 0.0
  AddAD1 = 0.0
  AddAD2 = 0.0
  AddAD3 = 0.0
  AddAD4 = 0.0
  Ovipos = 0.0
  NOvipos = 0.0
  
  ### compute mortalities ###
  
  # water- and density-dependent juvenile mortality
  ddmortLP = calc_mosq_atr(import3_Wet, parms$LP_Mortality)
  LP = LP/(1.0 + ddmortLP*sum(LP)*da*time)
  
  ## Kill juveniles if intervention day
  if(length(parms$LP_intervention_days) > 0){
    intervention_effect_LP = calc_intervention_percent(daystep, intervention_function_type = parms$LP_intervention_function, 
                                                       application_days_arg = as.numeric(parms$LP_intervention_days),
                                                       percent_effective_arg = as.numeric(parms$LP_intervention_percent),
                                                       max_duration_arg = as.numeric(parms$LP_intervention_max_duration),
                                                       decay_exponent_arg = as.numeric(parms$LP_intervention_decay))
    
    LP = LP * (1 - intervention_effect_LP)
    
    if(daystep %in% as.numeric(parms$LP_intervention_days)){
      # browser()
    }
  }
  
  # temperature-dependent adult mortality
  if(is.null(parms$Adult_Mortality)){
    dimortAD = 1 / calc_mosq_atr(import1_Temp, parms$Lifespan)
  } else{
    dimortAD = calc_mosq_atr(import1_Temp, parms$Adult_Mortality)
  }
  
  
  AD1 = AD1*exp(-dimortAD*time)
  AD2 = AD2*exp(-dimortAD*time)
  AD3 = AD3*exp(-dimortAD*time)
  AD4 = AD4*exp(-dimortAD*time)
  
  ## Kill adults if intervention day
  if(length(parms$Adult_intervention_days) > 0){
    intervention_effect_adult = calc_intervention_percent(daystep, intervention_function_type = parms$Adult_intervention_function, 
                                                          application_days_arg = as.numeric(parms$Adult_intervention_days),
                                                          percent_effective_arg = as.numeric(parms$Adult_intervention_percent),
                                                          max_duration_arg = as.numeric(parms$Adult_intervention_max_duration),
                                                          decay_exponent_arg = as.numeric(parms$Adult_intervention_decay))
    
    if(daystep %in% as.numeric(parms$Adult_intervention_days)){
      AD1 = AD1 * (1 - intervention_effect_adult)
      AD2 = AD2 * (1 - intervention_effect_adult)
      AD3 = AD3 * (1 - intervention_effect_adult)
      AD4 = AD4 * (1 - intervention_effect_adult)
    }
  }

  ### compute development rates ###
  
  # proportion of adults in diapause
  # DiapProp = calc_mosq_atr(import2_DayHrs, parms$Diapause)
  DiapProp = 0.0;
  
  # temperature-dependent egg development rate
  # if(import1_Temp > as.numeric(parms$Cutoff)){
  
  if(import1_Temp > as.numeric(parms$Cutoff) & import3_Wet > parms$Egg_cutoff){
    vtE = calc_mosq_atr(import1_Temp, parms$Egg)
  } else {
    vtE = 0.0
  }
  
  # temperature-dependent juvenile development rate
  # if(import1_Temp > as.numeric(parms$Cutoff) & import1_Temp <= as.numeric(parms$LP$Tm)){
  if(import1_Temp > as.numeric(parms$Cutoff)){
    vtLP = calc_mosq_atr(import1_Temp, parms$LP)
  } else {
    vtLP = 0.0
  }
  
  # temperature-dependent adult development rate
  if(import1_Temp > as.numeric(parms$Cutoff)){
    if(is.null(parms$GP)){
      vtA = 1 / calc_mosq_atr(import1_Temp, parms$Lifespan)
    } else{
      vtA = calc_mosq_atr(import1_Temp, parms$GP)
    }
  } else {
    vtA = 0.0
  }
  
  # reset the number of individuals that have reached each gonotrophic cycle to zero
  g1 = 0.0; g2 = 0.0; g3 = 0.0; g4 = 0.0
  
  ### shifting of age (development) bins ###
  if(import1_Temp > as.numeric(parms$Cutoff)){
    
    if(import3_Wet < parms$Egg_cutoff & parms$Egg_velocity_distribution != "constant"){
      MED = rep(0,length(ageEggs))  
      
      NewEggs = ConvolveFunc(x1 = Eggs, y1 = MED, padsize = length(Eggs))
    } else if(parms$Egg_velocity_distribution == "exponential"){
      MED = dexp(x = ageEggs, rate = 1/vtE*time)

      MED = MED/sum(MED)
      
      NewEggs = ConvolveFunc(x1 = Eggs, y1 = MED, padsize = length(Eggs))
    } else if(parms$Egg_velocity_distribution == "gamma"){
      MED = dgamma(x = ageEggs, shape = vtE*time, rate = 1)
      
      MED = MED/sum(MED)
      
      NewEggs = ConvolveFunc(x1 = Eggs, y1 = MED, padsize = length(Eggs))
    } else if(parms$Egg_velocity_distribution == "constant"){
      EggsShifted = ageEggs + rep(vtE*time, length(ageEggs))
    } else{
      print("Undefined egg velocity distribution")
      return(NA)
    }
    
    # if(sum(NewEggs$conv < 0) > 0 | NewEggs$emigrants < 0){
    #   browser()
    # }
    
    LPShifted = ageLP + rep(vtLP*time, length(ageLP))
    
    # Pupae that exceed the pupal break point enter the adult stage as below. All others
    # are shifted to the right with boxes replaced on the left side of the age (development) domain
    LP = c(rep(0.0, length(LPShifted[LPShifted > 1.0])), LP[LPShifted <= 1.0])
    
    # Eggs that exceed the egg breakpoint are added to the LP stage
    # To prevent spikes, we divide the new LP individuals evenly among the new bins
    if(length(parms$LP_intervention_days) > 0 | as.numeric(parms$LP_intervention_percent) > 0){
      if(daystep %in% as.numeric(parms$LP_intervention_days)){
        LP = LP * (1 - as.numeric(parms$LP_intervention_percent))
      }
    }
    
    if(parms$Egg_velocity_distribution == "exponential" | parms$Egg_velocity_distribution == "gamma"){
      lp1 = NewEggs$emigrants
      
      if(length(parms$LP_intervention_days) > 0){
        if(daystep %in% as.numeric(parms$LP_intervention_days)){
          lp1 = lp1 * (1 - as.numeric(parms$LP_intervention_percent))
        }
      }
      
      LP[1:length(LPShifted[LPShifted > 1.0])] = LP[1:length(LPShifted[LPShifted > 1.0])] +
        lp1/length(LPShifted[LPShifted > 1.0])
      
      AddLP = lp1
      
    } else if(parms$Egg_velocity_distribution == "constant"){
      lp1 = sum(Eggs[EggsShifted > 1.0])
      
      if(length(parms$LP_intervention_days) > 0){
        if(daystep %in% as.numeric(parms$LP_intervention_days)){
          lp1 = lp1 * (1 - as.numeric(parms$LP_intervention_percent))
        }
      }
      
      LP[1:length(LPShifted[LPShifted > 1.0])] = LP[1:length(LPShifted[LPShifted > 1.0])] +
        lp1/length(LPShifted[LPShifted > 1.0])
      
      AddLP = lp1
      
    } else{
      print("Undefined egg velocity distribution")
      return(NA)
    }
    
    # First gonotrophic cycle:
    # Computing the exponentially distributed development rate
    # this was formerly gamma distributed 
    # Exponential distribution is a type of gamma distribution
    
    # MAD = dgamma(x = ageAD1, shape = vtA*time, rate = 1)
    MAD = dexp(x = ageAD1, rate = 1/vtA*time)   # density of exponential distribution, not random
    MAD = MAD/sum(MAD) # Normalizing
    
    NewAD1 = ConvolveFunc(x1 = AD1, y1 = MAD, padsize = length(AD1))
    AD1 = NewAD1$conv
    
    ad1 = sum(LP[LPShifted > 1.0])
    
    ## Kill new adults if intervention day
    if(length(parms$Adult_intervention_days) > 0){
      if(daystep %in% as.numeric(parms$Adult_intervention_days)){
        ad1 = ad1 * (1 - as.numeric(parms$Adult_intervention_percent))
      }
    }

    AD1[1] =  AD1[1] + ad1 # dont divide by da here because we also multiply by da
    g1 = NewAD1$emigrants*da # multiply by da because this is the Reimann sum approximation
    
    AddAD1 = ad1
    
    # second gonotrophic cycle:
    NewAD2 = ConvolveFunc(x1 = AD2, y1 = MAD, padsize = length(AD2))
    AD2 = NewAD2$conv
    
    ad2 = g1/da
    
    AD2[1] =  AD2[1] + g1/da # divide by da because this is a Dirac delta function approximation
    g2 = NewAD2$emigrants*da # multiply by da because this is the Reimann sum approximation
    
    AddAD2 = ad2
    
    # Third gonotrophic cycle:
    NewAD3 = ConvolveFunc(x1 = AD3, y1 = MAD, padsize = length(AD3))
    AD3 = NewAD3$conv
    
    ad3 = g2/da
    
    AD3[1] =  AD3[1] + ad3 # divide by da because this is a Dirac delta function approximation
    g3 = NewAD3$emigrants*da # multiply by da because this is the Reimann sum approximation
    
    AddAD3 = ad3
    
    # fourth gonotrophic cycle:
    NewAD4 = ConvolveFunc(x1 = AD4, y1 = MAD, padsize = length(AD4))
    AD4 = NewAD4$conv
    
    ad4 = g3/da
    
    AD4[1] =  AD4[1] + ad4 # divide by da because this is a Dirac delta function approximation
    g4 = NewAD4$emigrants*da # multiply by da because this is the Reimann sum approximation
    
    AddAD4 = ad4

    # Eggs that exceed the break point enter the larval stage as above. All others are
    # shifted to the right with boxes replaced on the left side of the age (development) domain
    
    # Eggs are added to the first egg age (development) bin by adult oviposition
    # Ovipos = as.numeric(parms$Oviposition)*(g1 + g2 + g3 + g4)*time
    
    Ovipos = as.numeric(calc_mosq_atr(import1_Temp,parms$Oviposition))*(g1 + g2 + g3 + g4)*time
    
    # if(daystep == 20){
    #   print("oviposition value")
    #   print(parms$Oviposition)
    # }
    
    NOvipos = (g1 + g2 + g3 + g4)
    
    ### Sigmoidal function between x and y
    # Ovipos = as.numeric(sigmoidal(calc_mosq_atr(import1_Temp,parms$Oviposition),0.1,6))*(g1 + g2 + g3 + g4)*time
        
    ### Exponentiate oviposition so it is always positive
    
    # Ovipos = as.numeric(exp(calc_mosq_atr(import1_Temp,parms$Oviposition)))*(g1 + g2 + g3 + g4)*time
    
    if(parms$Egg_velocity_distribution == "exponential" | parms$Egg_velocity_distribution == "gamma"){
      Eggs = NewEggs$conv
      
      Eggs = Eggs + Ovipos/(length(Eggs)*da)
      
      AddEggs = (Ovipos/(length(Eggs)*da))*length(Eggs)
      
    } else if(parms$Egg_velocity_distribution == "constant"){
      Eggs = c(rep(0.0, length(EggsShifted[EggsShifted > 1.0])), Eggs[EggsShifted <= 1.0])
      
      Eggs[1:length(EggsShifted[EggsShifted > 1.0])] = Eggs[1:length(EggsShifted[EggsShifted > 1.0])] +
          Ovipos/(length(EggsShifted[EggsShifted > 1.0])*da)
      
      AddEggs = (Ovipos/(length(EggsShifted[EggsShifted > 1.0])*da))*length(EggsShifted[EggsShifted > 1.0])
    } else{
      print("Undefined egg velocity distribution")
      return(NA)
    }
    
  }
  
  # estimate of the number of active mosquitos
  ActivMosq = (1.0 - DiapProp)*(g1 + g2 + g3 + g4)
  TotalMosq = sum(AD1) + sum(AD2) + sum(AD3) + sum(AD4)
  
  # update state variables
  out = list(yinout = c(Eggs, LP, AD1, AD2, AD3, AD4), 
             data = c('ActMosq' = ActivMosq, 'TotMosq' = TotalMosq, 
                      'Eggs' = sum(Eggs), 'LP' = sum(LP),
                      'AD1' = sum(AD1), 'AD2' = sum(AD2), 'AD3' = sum(AD3), 'AD4' = sum(AD4),
                      'AddEggs' = AddEggs, 'AddLP' = AddLP, 'AddAD1' = AddAD1, 'AddAD2' = AddAD2, 'AddAD3' = AddAD3, 'AddAD4' = AddAD4,
                      'MortLP' = 1/(1.0 + ddmortLP*sum(LP)*da*time), 'MortAD' = exp(-dimortAD*time),
                      'Ovipos' = Ovipos, 'NOvipos' = NOvipos,
                      'vtE' = vtE, 'vtLP' = vtLP, 'vtA' = vtA))
  return(out)
}