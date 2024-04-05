[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10931935.svg)](https://doi.org/10.5281/zenodo.10931935)

# FeatureSelectionWithGammaMetric
Repository with the R scripts used for the simulation study of the feature selection with the gamma-metric.

# File functionGammaMetric.R
R scripts with all the functions needed to compute the $\gamma$-metric value of a set of features. The main function to call in other scripts is gammaMetric().

# File functionGeneration.R
This file contains the functions to generate data for all the scenarios. With the following description
- generation(): to generate the data of scenario 1 and 2.
- genData(): to generate the data of scenario 3
Two other functions, produceRho() and produceBeta() are used to generate the variance-covariance matrix (Sigma) and beta vector of scenario 3

# File functions.R
Contains a list of functions to call to apply each feature selection method (one for each method, one for BASELINE and one for BEST)
The function performSimulationIteration() which is called inside the parallel loop of the simulation. This function, when called, execute the feature selection (on train), 
apply the logistic regression  model (on train), compute the prediction of the model and the performance indicators for a given method. 
The dataframeSelectedFeaturesScenario2() and dataframeSelectedFeaturesScenario3() are two functions used to aggregate the results in a dataframe for the plots of the manuscrit.

# MainScenario1
The main script for scenario 1 in the simulation study. With definition of the generation parameters, parallel looping, saving results, aggregations of the results and plots.

# MainScenario2
The main script for scenario 2 in the simulation study. With definition of the generation parameters, parallel looping, saving results, aggregation of the results and plots.

# MainScenario3 
The main script for scenario 3 in the simulation study. With definition of the generation parameters, parallel looping, saving and aggregation of the results and plots.
