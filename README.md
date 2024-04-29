[![DOI](https://zenodo.org/badge/782467101.svg)](https://zenodo.org/doi/10.5281/zenodo.10931934)

# FeatureSelectionWithGammaMetric
Repository with the R scripts used for the simulation study of the feature selection with the $\gamma$-metric. The code was run on [SYSTEM SPEC]. Each scenario (MAinScenario1, MainScenario2, and MainScenario3) can be run independently of any other scenario. If one wish to execute one of the MainScenario make sure that all the files and folders of this repository are installed in the same directory and to set it as the working directory of your R environment. 

# Folder Scenario 1
Folder with all the results of scenario 1 of the simulation study (1 file).

# Folder Scenario 2
Folder with all the results of scenario 2 of the simulation study (16 files, one for each combination balance/separation).

# Folder Scenario 3 
Folder with all the results of scenario 3 of the simulation study (6 files, one for each constant/non-constant correlation with different correlation levels).

# File functionGammaMetric.R
R scripts with all the functions needed to compute the $\gamma$-metric value of a set of features. The main function to call to compute the $\gamma$-metric value of a set of features is gammaMetric().

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

# File MainScenario1.R
The main script for scenario 1 in the simulation study. With definition of the generation parameters, parallel looping, saving results, aggregations of the results and plots.
In the scrit, results files are saved in a folder called 'Scenario1' and the figures in a folder called 'Figures'. When running the script, saving are integrated in the loops. It is necessary then to change the path where the files are saved if one wish to generate new results without erasing the previous results. [TIME TO RUN WITH SYSTEM SPEC]. To run this script on less repetitions one can change the parameter called R in the script.

# File MainScenario2.R
The main script for scenario 2 in the simulation study. With definition of the generation parameters, parallel looping, saving results, aggregation of the results and plots.
In the scrit, results files are saved in a folder called 'Scenario2' and the figures in a folder called 'Figures'. When running the script, saving are integrated in the loops. It is necessary then to change the path where the files are saved if one wish to generate new results without erasing the previous results. [TIME TO RUN WITH SYSTEM SPEC]. To run this script on less repetitions one can change the parameter called R in the script.

# File MainScenario3.R 
The main script for scenario 3 in the simulation study. With definition of the generation parameters, parallel looping, saving and aggregation of the results and plots.
In the scrit, results files are saved in a folder called 'Scenario3' and the figures in a folder called 'Figures'. When running the script, saving are integrated in the loops. It is necessary then to change the path where the files are saved if one wish to generate new results without erasing the previous results. [TIME TO RUN WITH SYSTEM SPEC]. To run this script on less repetitions one can change the parameter called R in the script.

# File appendixFigures.R
This script contains the R code necessary to generate the figures and tables of the appendix and the supporting information. 
