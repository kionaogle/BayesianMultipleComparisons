# Bayesian Multiple Comparisons

A project aimed at identifying the need to account for multiple comparisons in hierarchical Bayesian models. The files contain code for running the four real-data examples in the main manuscript, and for calculating PPI (partial pooling index) based on the output from the JAGS simulations.

**Author names: Kiona Ogle, Drew Peltier, Michael Fell, Jessica Guo, Heather Kropp, Jarrett Barber**

Supplemental materials for Ogle et al. *(in prep)*.


# Explanation of folder structure and associated files

The **RealExamples** folder contains subfolders for each of the four real-data examples (*AllometricScaling*, *DogLearning*, *OrangeTrees*, and *PlantWaterStress*).

A master R script file is located within each example folder (e.g., *AllometricScaling_script.R*, *DogLearning_script.R*, etc.). These script files read-in data, create data lists for JAGS, set-up of the model variants (HB, CP, and NH) to run in JAGS, run the JAGS models, and compute PPI based on MCMCM output from each of the model variants. Users should begin by opening and browsing these script files.

The master R script files call model files that have been set-up for each of the model variants (HB, CP, and NH). Users can find the model files in the **models** subfolder that is contained within each of the example subfolders (e.g., *.../DogLearning/models*).

Each example sufolder also contains a subfolder with the data (*.../data*) and a subfolder (*.../source*) with the source code to compute PPI, which is called from the master R script; some examples subfolders may also contain a subfolder with initials used to initialize the JAGS models (e.g., *.../inits*), and all example subfolders contain an output (*.../output*) subfolder for storing results from the JAGS simulations upon executing the master R scripts (as posted, the output subfolders do not contain data, just a README file).


[![DOI](https://zenodo.org/badge/153157527.svg)](https://zenodo.org/badge/latestdoi/153157527)
