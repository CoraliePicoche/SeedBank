This folder contains all scripts used to model the effect of seed banks on phytoplanktonic communities.  

* `beta_law.r` only plots a figure for a potential distribution of sinking rate
* `diagnostics_single_simu.r` plots different representations of a simulation: variation in parameters due to quadratic programming, time series of simulated abundances, comparison of simulated growth rates to observed growth rates
* `generate_fixed_values.r` is as simple generator of species-specific parameter values for sinking rates and niche widths, to remove this source of variability in the rest of the simulations
* `infer_interaction_matrix_growth_rate.r` contains the functions to turn MAR coefficients to BH coefficients, and to calibrate them with quadratic programming
* `matrix_MAR_clean.r` is a function that takes a MARSS object (interaction matrix determined by a multivariate autoregressive model in Picoche & Barraquand 2020) into an easier-to-treat matrix
* `phenology.r` uses observed time series to infer the status (generalist/specialist) and thermal optimum of each species. It also plots annual blooms and statistics on their number, duration, and temperature trigger
* `step_functions.r` contains all the functions for each time step, i.e. growth step and exchange step
* `summary_statistics.r` computes the statistics used in the calibration of the interactions
* `thermal_niche.r` draws the growth response of analysed species to temperature, depending on their niche width and thermal optimum
* `working_on_temperature.r` compares modeled and observed temperatures

Scripts for sensitivity analysis of the models and test of different scenarios without seed banks can be found in `common_scripts.r`.
