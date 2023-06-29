# Stability and performance evaluation of different windspeed spatial interpolation models to assist pesticide dispersion estimates.  

## Contents

```
└───thesis_project_f2
    │   README.md
    │   [1] spatial-interpolation-overview.Rmd
    │   [2] dynamic-training-routine.R
    │   [3] mqrbf-model-simple.R
    │   [4] thesis-validation-overview.Rmd
    └───Data
    	│   [ 5] weights_with_lowest_rmse_per_opt.rds
    	│   [ 6] knmi_windspeed_complete_set.rds
    	│   [ 7] knmi_windspeed_observations_training_set.rds
    	│   [ 8] knmi_windspeed_receptors_training_set.rds
    	│   [ 9] dynamic_training_routine_training_results_annual_avg.rds
    	│   [10] dynamic_training_routine_training_settings_annual_avg.rds
    	│   [11] source_map_crs.rds
    	│   [12] netherlands_reference_map.gpkg
    	│   [13] knmi_weather_stations.csv
    	│   [14] knmi_winspeed_days_2017.json
	
    └───Docs 
    	│   [15] dynamic-training-routine.md
    	│   [16] Thesis paper.pdf
    	│   [17] Thesis paper.pdf
    	│   [18] Spatial statistics overview.pdf
    	│   [19] spatial-interpolation-overview.html
    	│   [20] thesis-validation-overview.html
```

## Description

[1]: Contains the R markdown notebook explaining the methods and thinking process used for the different spatial interpolation methods. This notebook contains a complete overview of the steps but not the code used for training of the weights.

[2]: Contains the R code for the dynamic training routine used to find ideal model hyper parameters for the spatial interpolation models used and storing the results. 

[3]: Contains a simple implementation of the MQ-RBF model. 

[4]: Contains the R markdown notebook explaining and showing the results of the different validation metrics used for the spatial interpolation methods using results from the dynamic training routine.

[5]: File containing BAP hyperparameters.  

[6]: File containing annual observations as SF object.

[7]: File containing observation of `01-01-2017` as SF object used to calculate BDP.

[8]: File containing receptors as SF object.

[9]: File containing BDP hyperparameters.

[10]: File containing the settings used to calculate BDP hyperparameters.

[11]: File containing the CRS of BRP Gewaspercelen map used as source map.

[12]: File containing the reference map used to perform tesselation.

[13]: File containing the coordinates in (WG:84) of the KNMI stations.

[14]: File containing daily weather data for windspeed of the KNMI stations in the year 2017. 

[15]: Contains the documentation how to use and tweak the dynamic training routine of file [2].

[16]: Copy of thesis paper

[17]: Copy of thesis poster used in presentation.

[18]: Copy of asside used in appendix of thesis paper.

[19]: Knitted html file of execution in file [1].

[20]: Knitted html file of execution in file [4]
