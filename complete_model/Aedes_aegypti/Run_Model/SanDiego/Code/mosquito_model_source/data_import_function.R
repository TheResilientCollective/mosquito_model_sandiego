data.import <- function(config_data, burnin = 14){
  
  # read in mosquito trap data
  observed.data <- read.csv(config_data$MOSQ_DATA_PATH)
  
  # read in weather data
  ExternalInputs = read.csv(config_data$EXT_DATA_PATH, h = T) 
  
  # browser()
  
  earliest_common_daystep = max(c(min(observed.data$Daystep), min(ExternalInputs$Daystep)))
  
  observed.data <- observed.data %>%
    dplyr::filter(Daystep > earliest_common_daystep)
  
  ExternalInputs <- ExternalInputs %>%
    dplyr::filter(Daystep > earliest_common_daystep)
  
  ### Training on a selected % of the seasons ###
  if(!is.null(config_data$TRAIN_PCT)){
    train_pct <- as.numeric(config_data$TRAIN_PCT)
  } else{
    train_pct <- 0.8
  }
  
  if(config_data$MODEL_RUN_TYPE == "train"){
    
    # calculate the number of seasons
    observed.data <- observed.data %>%
      mutate(Season = trunc(observed.data$Daystep/365)) %>%
      mutate(Set = ifelse(Season %in% unique(Season)[1:trunc(train_pct*length(unique(Season)))], 'train', 'test')) %>%
      #mutate(Set = ifelse(Season %in% unique(Season)[1:trunc(.50*length(unique(Season)))], 'train', 'test')) %>%
      dplyr::filter(Set == 'train') %>%
      dplyr::select(-Set, -Season)
  }
  
  # burn-in period included in both training and testing sets
  daystep_start <- min(observed.data$Daystep) - burnin
  daystep_end <- ifelse(config_data$MODEL_RUN_TYPE == "train", max(observed.data$Daystep), max(ExternalInputs$Daystep)) 
  
  daystep_vec <- daystep_start:daystep_end
  
  # rename 'Daystep' column to 'Datestep'
  ExternalInputs <- ExternalInputs %>%
    dplyr::rename(Datestep = Daystep)
  
  ### handle water data ###
  wet_cols <- config_data$WETNESS$column  # water data variables
  wet_weights <- config_data$WETNESS$val  # weights of water data variables
  
  # weighted water input
  wet_list <- list()
  for (i in 1:length(wet_cols)){
    wet_list[i] <- ExternalInputs[wet_cols[i]] * wet_weights[i]
  }
  
  ExternalInputs$Wetness <- Reduce('+', wet_list)
    
  external.data <- ExternalInputs %>%
    dplyr::filter(Datestep %in% daystep_vec) %>%
    
    # normalize daylight hours
    mutate(DaylightHrs = Daylighthrs - max(Daylighthrs)) %>%
    mutate(TimeSeq = 0:(max(Datestep) - min(Datestep))) %>%
    dplyr::rename(Daystep = Datestep, Temperature = TempC) %>%
    dplyr::select(-Daylighthrs)
  
  data_list <- list('observed' = observed.data, 'external'= external.data)
  return(data_list) 
}
