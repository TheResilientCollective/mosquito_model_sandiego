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
  if(import1_Temp < 0.0){
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
  
  ### compute mortalities ###
  
  # water- and density-dependent juvenile mortality
  ddmortLP = calc_mosq_atr(import3_Wet, parms$LP_Mortality)
  LP = LP/(1.0 + ddmortLP*sum(LP)*da*time)
  
  # temperature-dependent adult mortality
  dimortAD = calc_mosq_atr(import1_Temp, parms$Adult_Mortality)
  AD1 = AD1*exp(-dimortAD*time)
  AD2 = AD2*exp(-dimortAD*time)
  AD3 = AD3*exp(-dimortAD*time)
  AD4 = AD4*exp(-dimortAD*time)

  ### compute development rates ###
  
  # proportion of adults in diapause
  DiapProp = calc_mosq_atr(import2_DayHrs, parms$Diapause)
  
  # temperature-dependent egg development rate
  if(import1_Temp > as.numeric(parms$Cutoff)){
    vtE = calc_mosq_atr(import1_Temp, parms$Egg)
  } else {
    vtE = 0.0
  }
  
  # temperature-dependent juvenile development rate
  if(import1_Temp > as.numeric(parms$Cutoff) & import1_Temp <= as.numeric(parms$LP$Tm)){
    vtLP = calc_mosq_atr(import1_Temp, parms$LP)
  } else {
    vtLP = 0.0
  }
  
  # temperature-dependent adult development rate
  if(import1_Temp > as.numeric(parms$Cutoff)){
    vtA = calc_mosq_atr(import1_Temp, parms$GP)
  } else {
    vtA = 0.0
  }
  
  # reset the number of individuals that have reached each gonotrophic cycle to zero
  g1 = 0.0; g2 = 0.0; g3 = 0.0; g4 = 0.0
  
  ### shifting of age (development) bins ###
  if(import1_Temp > as.numeric(parms$Cutoff)){
    EggsShifted = ageEggs + rep(vtE*time, length(ageEggs))
    LPShifted = ageLP + rep(vtLP*time, length(ageLP))
    
    # Pupae that exceed the pupal break point enter the adult stage as below. All others
    # are shifted to the right with boxes replaced on the left side of the age (development) domain
    LP = c(rep(0.0, length(LPShifted[LPShifted > 1.0])), LP[LPShifted <= 1.0])
    
    # Eggs that exceed the egg breakpoint are added to the LP stage
    # To prevent spikes, we divide the new LP individuals evenly among the new bins
    LP[1:length(LPShifted[LPShifted > 1.0])] = LP[1:length(LPShifted[LPShifted > 1.0])] +
      sum(Eggs[EggsShifted > 1.0])/length(LPShifted[LPShifted > 1.0])
    
    # First gonotrophic cycle:
    # Computing the exponentially distributed development rate
    # this was formerly gamma distributed 
    # Exponential distribution is a type of gamma distribution
    
    # MAD = dgamma(x = ageAD1, shape = vtA*time, rate = 1)
    MAD = dexp(x = ageAD1, rate = 1/vtA*time)   # density of exponential distribution, not random
    MAD = MAD/sum(MAD) # Normalizing
    
    NewAD1 = ConvolveFunc(x1 = AD1, y1 = MAD, padsize = length(AD1))
    AD1 = NewAD1$conv
    AD1[1] =  AD1[1] + sum(LP[LPShifted > 1.0]) # dont divide by da here because we also multiply by da
    g1 = NewAD1$emigrants*da # multiply by da because this is the Reimann sum approximation
    
    # second gonotrophic cycle:
    NewAD2 = ConvolveFunc(x1 = AD2, y1 = MAD, padsize = length(AD2))
    AD2 = NewAD2$conv
    AD2[1] =  AD2[1] + g1/da # divide by da because this is a Dirac delta function approximation
    g2 = NewAD2$emigrants*da # multiply by da because this is the Reimann sum approximation
    
    # Third gonotrophic cycle:
    NewAD3 = ConvolveFunc(x1 = AD3, y1 = MAD, padsize = length(AD3))
    AD3 = NewAD3$conv
    AD3[1] =  AD3[1] + g2/da # divide by da because this is a Dirac delta function approximation
    g3 = NewAD3$emigrants*da # multiply by da because this is the Reimann sum approximation
    
    # fourth gonotrophic cycle:
    NewAD4 = ConvolveFunc(x1 = AD4, y1 = MAD, padsize = length(AD4))
    AD4 = NewAD4$conv
    AD4[1] =  AD4[1] + g3/da # divide by da because this is a Dirac delta function approximation
    g4 = NewAD4$emigrants*da # multiply by da because this is the Reimann sum approximation
    
    # Eggs that exceed the break point enter the larval stage as above. All others are
    # shifted to the right with boxes replaced on the left side of the age (development) domain
    Eggs = c(rep(0.0, length(EggsShifted[EggsShifted > 1.0])), Eggs[EggsShifted <= 1.0])
    
    # Eggs are added to the first egg age (development) bin by adult oviposition
    Ovipos = as.numeric(parms$Oviposition)*(1.0 - DiapProp)*(g1 + g2 + g3 + g4)*time

    Eggs[1:length(EggsShifted[EggsShifted > 1.0])] = Eggs[1:length(EggsShifted[EggsShifted > 1.0])] + 
      Ovipos/(length(EggsShifted[EggsShifted > 1.0])*da)
  }
  
  # estimate of the number of active mosquitos
  ActivMosq = (1.0 - DiapProp)*(g1 + g2 + g3 + g4)
  TotalMosq = sum(AD1) + sum(AD2) + sum(AD3) + sum(AD4)
  
  # update state variables
  out = list(yinout = c(Eggs, LP, AD1, AD2, AD3, AD4), 
             data = c('ActMosq' = ActivMosq, 'TotMosq' = TotalMosq, 
                      'Eggs' = sum(Eggs), "LP" = sum(LP)))
  return(out)
}