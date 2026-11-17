# Running fitted mosquito model

The scripts in this directory run the fitted mosquito model. 

Scripts:

- run_model.R: Runs one iteration of the fitted model and outputs time series of modeled population
    - `Rscript run_model_1rep.R alpha1_value beta1_value life_stage config_file_path`
      - alpha1_value and beta1_value are decimal numbers
      - life_stage is one of 4 options: "ActMosq" (active, biting), "TotMosq" (adults, non biting), "LP" (juveniles), "Eggs"
    - Inputs:
      - Put in the ./Input/ directory
        - ```fitted_params.csv```: this is FIT_PARAM_DATA_PATH in ```simulator_config.yaml```. This is a .csv file containing the names and fits of the fitted juvenile mortality parameters and initial conditions. NOTE: this does not need to be altered/edited and is used to format the inputs to the mosquito simulator. The real values used for alpha1 and beta1 are given as function arguments to the mosq_sim() function
        - Observed mosquito data: this is the MOSQ_DATA_PATH in ```simulator_config.yaml```. Contains columns "Daystep", "muMosq", "cost_weight" (this weights the target function for fitting and is set to 0.1 in the mosquito 'off season' and 1 in the 'on season'), and "Date". The columns must have those names.
        - Weather data: this is EXT_DATA_PATH in ```simulator_config.yaml```. Must have the columns "Daystep", "TempC", and "Date". The remaining columns are weather data, averaged for each day. The names in the weather columns must match the WETNESS: column: arg in ```simulator_config.yaml```
- plot_model_intervals.R: Samples the posterior distributions of alpha1 and beta1 (juvenile mortality parameters) to generate intervals for estimated mosquito populations