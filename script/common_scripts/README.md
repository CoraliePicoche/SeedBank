This folder contains used for the main analyses of the manuscript. 

* `compet_growth_rates_several_conditions.r` outputs the realized growth rate of species depending on the temperature and abundances of the other species (Section S7 of the SI).
* `main_sensitivity.r` runs the sensitivity analysis (including Fig. 3 in manuscript).
* `no_seed_bank_compet.r` launches simulations with and without a seed bank when interaction strengths are varied (Fig. 4 in manuscript) and compares the probabilities of survival of different species as a function of their characteristics (Fig. 5 in manuscript).  In `no_seed_bank_compet_keeping_intra.r`, the intraspecific interaction strengths are not modified while interspecific interactions strengths are (Section S6 in SI).
* `no_seed_bank_temp.r` launches simulations with and without a seed bank varying mean and variance of the temperature (Fig. 6 in manuscript).
* `step_functions.r` contains the main functions for temperature-dependent growth rate and the two substeps of the model (growth and exchanges).
