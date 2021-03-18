This folder contains all scripts used to model the effect of seed banks on phytoplanktonic communities.  

* `beta_law.r` only plots a figure for a potential distribution of sinking rate (explained in Section S3 in SI).
* `diagnostics_single_simu.r` plots different representations of a simulation: time series of simulated abundances (Fig. 2 in manuscript + Fig. S5 in SI), comparison of simulated growth rates to observed growth rates.
* `generate_fixed_values.r` is as simple generator of species-specific parameter values for sinking rates and niche widths, to remove this source of variability in the rest of the simulations.
* `growth_rate_presentation.r` shows the general functioning of the temperature-dependent growth rate. This is shown in Section S1 of SI.
* `infer_interaction_matrix_growth_rate.r` contains the function to turn MAR coefficients to BH coefficients.
* `matrix_MAR_clean.r` is a function that takes a MARSS object (interaction matrix determined by a multivariate autoregressive model in Picoche & Barraquand 2020) into an easier-to-treat matrix.
* `phenology.r` uses observed time series to infer the status (generalist/specialist) and thermal optimum of each species. It also plots annual blooms and statistics on their number, duration, and temperature trigger (main results are explained in Section S1 in SI).
* `summary_statistics.r` computes the statistics used in the calibration of the interactions.
* `thermal_niche.r` draws the growth response of analysed species to temperature, depending on their niche width and thermal optimum (shown in Section S1 in SI).
* `working_on_temperature.r` compares modeled and observed temperatures.

Scripts for sensitivity analysis of the models and test of different scenarios without seed banks can be found in `common_scripts.r`.
