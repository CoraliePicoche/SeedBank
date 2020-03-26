This folder contains all scripts used to model the effect of seed banks on phytoplanktonic communities. It is 

#### Main
* `beta_law.r` only plots a figure for a potential distribution of sinking rate
* `comparison_growth_rate.r` plots different values and models from the literature for phytoplanktonic growth rates. It is completed by `variation_in_area.r`, which focuses on the impact of the niche sizeand the metabolism for a specific formula (Scranton & Vasseur 2016)
* `diagnostics_multiple_simu.r` compares several simulations, each described by a specific file, mainly based on a visual comparison simulated abudances (NOTE: it is not up-to-date right now 26/03/2020)
* `diagnostics_single_simu.r` plots different representations of a simulation: variation in parameters due to quadratic programming, time series of simulated abundances, comparison of simulated growth rates to observed growth rates
* `example_expected_results.r` only plots mock figures to represent our hypotheses and the way results could be presented
* `example_niches.r` plots simplistic scenarios of niche construction (instead of a mechanistic model, a simple polynome determines the niche)
* `generate_fixed_values.r` is as simple generator of species-specific parameter values for sinking rates and niche widths, to remove this source of variability in the rest of the simulations
* `infer_interaction_matrix_growth_rate.r` contains the functions to turn MAR coefficients to BH coefficients, and to calibrate them with quadratic programming
* `main.r` is the script that retrieves all parameter values, simulation file, and launches the simulation to output text files of the simulated time series
* `matrix_MAR_clean.r` is a function that takes a MARSS object (interaction matrix determined by a multivariate autoregressive model in Picoche & Barraquand 2020) into an easier-to-treat matrix
* `phenology.r` uses observed time series to infer the status (generalist/specialist) and thermal optimum of each species. It also plots annual blooms and statistics on their number, duration, and temperature trigger
* `step_functions.r` contains all the functions for each time step, i.e. growth step and exchange step

#### Exploratory analyses
