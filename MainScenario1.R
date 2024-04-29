###########################################################################################################
#                                                 SCENARIO 1                                              #
###########################################################################################################

## Environment --------------------------------------------------------------------------------------------
if(!require(devtools)){
  # Package used to install packages
  install.packages("devtools")
}

if(!require(FSelector)){
  # Package to call feature selection methods
  devtools::install_version('FSelector', version = '0.34')
  library(FSelector)
}

if(!require(caret)){
  # Package used to calibrate models (and compute some indicators)
  devtools::install_version('caret', version = '6.0.94')
  library(caret)
}

if(!require(pROC)){
  # Package used to compute performance indicators
  devtools::install_version('pROC', version = '1.18.4')
  library(pROC)
}

if(!require(stringr)){
  # Package used to manipulate strings
  devtools::install_version('stringr', version = '1.5.0')
  library(stringr)
}

if(!require(ggplot2)){
  # Package used to create figures
  devtools::install_version('ggplot2', version = '3.4.3')
  library(ggplot2)
}

if(!require(mltools)){
  # Package used to compute MCC
  devtools::install_version('mltools', version = '0.3.5')
  library(mltools)
}

if(!require(parallel)){
  # Pacakage used for parallel work
  devtools::install_version('parallel', version = '4.3.1')
  library(parallel)
}

if(!require(doSNOW)){
  # Package used for parallel work
  devtools::install_version('doSNOW', version = '1.0.20')
  library(doSNOW)
}

if(!require(doRNG)){
  # Package used for reproductibility of parallel executions (seed)
  devtools::install_version('doRNG', version = '1.8.6')
  library(doRNG)
}

if(!require(ggpubr)){
  # Package used for grids in plots
  devtools::install_version('ggpubr', version = '0.6.0')
  library(ggpubr)
}
#### Sources --------------------------------------------------------------------------------------------
source('functionGammaMetric.R')
source('functionGeneration.R')
source('functions.R')

#### Parameters of the generation -----------------------------------------------------------------------
N <- 2000                                 # Number of observations
n_var_inf <- 3                            # Number of informative features
n_noise <- 22                             # Number of non-informative features
beta <- c(3, -2, 0.5, rep(0, n_noise))    # Beta coefficients for the generation process
beta0 <- 0                                # Intercept value
R <- 50                                   # Number of repetitions
set.seed(123)                             # Seed of the simulation


#### Parameters of the parallel computations ------------------------------------------------------------

# Feature selection methods 
Approach <- c('BASELINE', 'BEST', 'CFS', 'CHI2', 'CONS', 'IG', 'IGR', 'ONER', 
              'RELIEF', 'RFI', 'SU', 'SVM-RFE', 'GAMMA_BACK', 'GAMMA_BF', 'GAMMA_FORW', 'GAMMA_HC')

# Number of clusters
n_clusters <- parallel::detectCores() - 1
cl <- parallel::makeCluster(n_clusters)
doSNOW::registerDoSNOW(cl)
time <- vector('numeric', R)
res <- NULL

# Progress bar
pb <- txtProgressBar(min = 0, max = R*length(Approach), style = 3)
progress <- function(n) setTxtProgressBar(pb, n + iteration_count)
opts <- list(progress = progress)
iteration_count <<- 0

#### BEGINING OF THE SIMULATION -------------------------------------------------------------------------

for(i in 1:R){
  # Time of computation
  t0 <- Sys.time()
  
  # Generation of the training data
  dat_train <- generation(n = N, n_var_inf = n_var_inf, n_noise = n_noise, beta = c(beta0, beta))
  X_train <- as.matrix(dat_train[, -1])
  Y_train <- as.matrix(dat_train[, 1])
  
  # Generation of the test data
  dat_test <- generation(n = N, n_var_inf = n_var_inf, n_noise = n_noise, beta = c(beta0, beta))
  X_test <- as.matrix(dat_test[, -1])
  Y_test <- as.matrix(dat_test[, 1])
  
  # Perform feature selection, modelisation and prediction in parallel for each feature selection methods
  resIteration <- foreach(fs = Approach, .combine = 'rbind', .export = ls()[!ls() %in% c('performSimulationIteration', 'X_train', 'Y_train', 'X_test', 'Y_test')], .options.snow = opts, .verbose = FALSE) %dorng% {
    localRes <- performSimulationIteration(X_train, Y_train, X_test, Y_test, beta = beta, fs)
    row.names(localRes) <- NULL
    return(localRes)
  }
  
  # Updates
  t1 <- Sys.time()
  iteration_count <<- i*length(Approach)
  res <- rbind(res, resIteration)
  time[i] <- as.numeric(difftime(t1, t0, units = 'secs'))
}
close(pb)
parallel::stopCluster(cl)

#### Saving the results ----------------------------------------------------------------------------------
file_path <- paste0('Scenario1/res_', R, '_iterations.txt')
write.table(res, file_path)

############################################################################################################################################################################
############################################################################################################################################################################

#### Uploading the results for graphics and tables -------------------------------------------------------
res <- read.table('Scenario1/res_50_iterations.txt')

# Formating
res <- data.frame(res)

# Aggregated results
RES <- data.frame(Approach = Approach)

for(i in 1:nrow(RES)){
  # Line
  lignes <- res[which(res$Approach == RES[i, 'Approach']), ]
  
  # Number of features
  RES[i, 'NbVarSelected'] <- mean(lignes[, 'NbVarSelected'])
  RES[i, 'NbVarInf'] <- mean(lignes[, 'NbVarInf'])
  
  # Selection indicators
  RES[i, 'Accuracy_selection'] <- mean(lignes[, 'Accuracy_selection'])*100
  RES[i, 'Specificity_selection'] <- mean(lignes[, 'Specificity_selection'])*100
  RES[i, 'Sensitivity_selection'] <- mean(lignes[, 'Sensitivity_selection'])*100
  RES[i, 'MCC_selection'] <- mean(lignes[, 'MCC_selection'])
  
  # Classification indicators (training data)
  RES[i, 'Accuracy_training'] <- mean(lignes[, 'Accuracy_training'])*100
  RES[i, 'Specificity_training'] <- mean(lignes[, 'Specificity_training'])*100
  RES[i, 'Sensitivity_training'] <- mean(lignes[, 'Sensitivity_training'])*100
  RES[i, 'MCC_training'] <- mean(lignes[, 'MCC_training'])
  
  # Classification indicators (test data)
  RES[i, 'Accuracy_test'] <- mean(lignes[, 'Accuracy_test'])*100
  RES[i, 'Specificity_test'] <- mean(lignes[, 'Specificity_test'])*100
  RES[i, 'Sensitivity_test'] <- mean(lignes[, 'Sensitivity_test'])*100
  RES[i, 'MCC_test'] <- mean(lignes[, 'MCC_test'])
}

## Code for table 2 --------------------------------------------------------------------------------------
print(RES[, c('Approach', "NbVarSelected", 'NbVarInf', "Specificity_selection", "Sensitivity_selection", "Accuracy_test")])

#### Plot of selected features ---------------------------------------------------------------------------

## Code for Figure 1 -------------------------------------------------------------------------------------

# Data frame with all the features selected by each method
VARSELECTED <- data.frame(sapply(Approach, FUN = function(x){
  approach_fs <- which(res$Approach == x)
  list_var_selected <- unlist(strsplit(res$VarSelected[approach_fs], ';'))
  indices <- match(c(paste0('x', 1:n_var_inf), paste0('Noise', 1:n_noise)), names(table(list_var_selected)))
  table(list_var_selected)[indices]
}))

# Names of the features as a column
VARSELECTED$Feature <- c(paste0('x', 1:n_var_inf), paste0('Noise', 1:n_noise))
colnames(VARSELECTED)[colnames(VARSELECTED) == 'SVM.RFE'] <- 'SVM-RFE' 

# Switch to a long format table
VARSELECTED_LONG <- reshape2::melt(VARSELECTED, id.vars = 'Feature')
plot_figure_01 <- ggplot(data = VARSELECTED_LONG, aes(x = Feature, y = variable))+
  geom_tile(aes(fill = value), color = 'white', lwd = 0.3, linetype = 1)+
  scale_x_discrete(breaks = c(paste0('x', 1:n_var_inf), paste0('Noise', 1:n_noise)), limits = c(paste0('x', 1:n_var_inf), paste0('Noise', 1:n_noise)), labels = c(paste0('x', 1:n_var_inf), paste0('N', 1:n_noise)))+
  scale_y_discrete(limits = rev(Approach))+
  scale_fill_gradientn(colours = rev(viridis::inferno(10)), breaks = seq(0, R, R/5))+
  labs(fill = '', x = '', y = '')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 30),
        axis.text.y = element_text(size= 35), 
        legend.text = element_text(size = 25),
        legend.key.size = unit(5, 'cm'), 
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(2, 'cm'))

setEPS()
postscript('Figures/Figure_01.eps', horizontal = FALSE, width = 18, height = 10)
plot_figure_01
dev.off()

## Code for Figure 2 -------------------------------------------------------------------------------------
plot_figure_02 <- ggplot(data = res, aes(x = Approach, y = Accuracy_test*100))+
  geom_boxplot(size = 0.3, outlier.size = 0.5, outlier.shape = 19) +
  coord_flip()+
  scale_x_discrete(limits = rev(Approach), labels = rev(Approach))+
  scale_y_continuous(limits = c(40, 100))+
  grids(axis = 'y', color = 'grey', linetype = 'dashed', size = 0.3)+
  grids(axis = 'x', color = 'grey', linetype = 'solid', size = 0.3)+
  labs(y = '', x = '')+
  theme(axis.text = element_text(size = 35),
        panel.background = element_rect(fill = 'white', colour = 'black', linewidth = 0.2))

setEPS()
postscript('Figures/Figure_02.eps', horizontal = FALSE, width = 18, height = 10)
plot_figure_02
dev.off()
