# Requirements
- Data: appropriately structured data on bicycle demand for CityBike (available on the course page, or from https://citibikenyc.com/system-data)
- R packages installed: c("ggplot2", "lubridate", "dplyr", "tidyr", "purrr", "tsibble", "fable", "feasts", "fabletools", "forecast","bootUR", "tidyverse", "FoReco")  

# Usage
- Four scripts are available
  - exploration.R replicates the data exploration
  - unit_root_and_ets.R replicates the unit-root testing and ETS model fitting and out-of-sample exercise
  - dummy_sarima.R replicates the dummy + sarima model fitting and out-of-sample exercise
  - forecast_recon.R replicates the forecast reconciliation and out-of-sample exercise

# Outputs
- plots in folder /output
- data for the out-of-sample comparison of the models in folder /output

# LICENSE 
MIT

