####################################################################################################################################
#                                                           SCENARIO 2                                                             #  
####################################################################################################################################

## Environment ---------------------------------------------------------------------------------------------------------------------
if(!require(devtools)){
  # Package used to install packages
  install.packages('devtools')
}

if(!require(FSelector)){
  # Package used for feature selection
  devtools::install_version('FSelector', version = '0.34')
  library(FSelector)
}

if(!require(caret)){
  # Package used to calibrate models (and compute indicators)
  devtools::package_version('caret', version = '6.0.94')
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

if(!require(parallel)){
  # Package to execute parallel work
  devtools::install_version('parallel', version = '4.3.1')
  library(parallel)
}

if(!require(doSNOW)){
  # Package used for parallelisation
  devtools::install_version('doSNOW', version = '1.0.20')
  library(doSNOW)
}

if(!require(doRNG)){
  # Package used for reproductibility of parallel looping (define a seed for parallel loop)
  devtools::install_version('doRNG', version = '1.8.6')
  library(doRNG)
}

if(!require(ggplot2)){
  # Package used for graphics
  devtools::install_version('ggplot2', version = '3.4.3')
  library(ggplot2)
}

if(!require(ggpubr)){
  # Package used for grids in ggplot graphs
  devtools::install_version('ggpubr', version = '0.6.0')
  library(ggpubr)
}


## Sources -------------------------------------------------------------------------------------------------------------------------
source('functionGammaMetric.R')
source('functionGeneration.R')
source('functions.R')

## Parameters of the generation ----------------------------------------------------------------------------------------------------
N <- 2000                                          # Number of observations
n_var_inf <- 3                                     # Number of informative features
n_noise <- c(47, 97, 147, 197)                     # Number of non-informative features
separation <- c('faible', 'forte')                 # Separation between classes
repartition <- c('equilibre', 'desequilibre')      # Repartition between classes
R <- 50                                             # Number of repetition
set.seed(123)                                      # Seed of the simulation

# Feature selection methods
Approach <- c('BASELINE', 'BEST', 'CFS', 'CHI2', 'CONS', 'IG', 'IGR', 'ONER', 
              'RELIEF', 'RFI', 'SU', 'SVM-RFE', 'GAMMA_BACK', 'GAMMA_BF', 'GAMMA_FORW', 'GAMMA_HC')

# EXECUTE TO RUN THE SIMULATIONS OR GO TO UPLOADING THE RESULTS
###########################################################################################################################################################################
###########################################################################################################################################################################
## Parameters of the parallel computations -----------------------------------------------------------------------------------------

# Number of clusters
n_clusters <- parallel::detectCores()-1
cl <- parallel::makeCluster(n_clusters)
doSNOW::registerDoSNOW(cl)

# Time vector
time <- expand.grid(1:R, repartition, separation, (n_noise + n_var_inf), 0)
colnames(time) <- c('Repetition', 'Repartition', 'Separation', 'NbVar', 'Time')

# Initialisations
iteration_count <<- 0
i <- 1
res <- NULL

# Progress bar
total_iteration <- R*length(Approach)*length(n_noise)*length(repartition)*length(separation)
pb <- txtProgressBar(min = 0, max = total_iteration, style = 3)
progress <- function(n) setTxtProgressBar(pb, n + iteration_count)
opts <- list(progress = progress)

## BEGINING OF THE SIMULATION -----------------------------------------------------------------------------------------------------
for(p in 1:length(n_noise)){
  for(sep in separation){
    for(repa in repartition){
      # Time starter
      t0 <- Sys.time()
      
      # Definition of the beta vector
      if(sep == 'faible'){
        # In the case of weak separation, balanced and unbalanced classes
        beta0 <- ifelse(repa == 'equilibre', 0.5, -2.75)
        beta <- c(0.6, -2.5, -1, rep(0, n_noise[p]))
      }else{
        if(repa == 'equilibre'){
          # In the case of strong separation and balanced classes
          beta0 <- 0
          beta <- c(3.6, -4, -1, rep(0, n_noise[p]))
        }else{
          # In the case of strong separation and unbalanced classes
          beta0 <- -2.65
          beta <- c(3.6, -2.2, -1, rep(0, n_noise[p]))
        }
      }
      
      for(r in 1:R){
        # Generation of the training data
        dat_train <- generation(n = N, n_var_inf = n_var_inf, n_noise = n_noise[p], beta = c(beta0, beta))
        X_train <- as.matrix(dat_train[, -1])
        Y_train <- as.matrix(dat_train[, 1])
        
        # Generation of the validation data
        dat_test <- generation(n = N, n_var_inf = n_var_inf, n_noise = n_noise[p], beta = c(beta0, beta))
        X_test <- as.matrix(dat_test[, -1])
        Y_test <- as.matrix(dat_test[, 1])
        
        # Perform feature selection, modelisation and prediction in parallel for each feature selection methods
        resIteration <- foreach(fs = Approach, .combine = 'rbind', .export = ls()[!ls() %in% c('performSimulationIteration', 'X_train', 'Y_train', 'X_test', 'Y_test', 'beta')], .options.snow = opts, .verbose = FALSE) %dorng% {
          localRes <- performSimulationIteration(X_train = X_train, Y_train = Y_train, X_test = X_test, Y_test = Y_test, beta = beta, fs = fs)
          row.names(localRes) <- NULL
          return(localRes)
        }
        
        # Filling the common column 
        resIteration$Separation <- sep
        resIteration$Repartition <- repa
        resIteration$NbVar = n_var_inf + n_noise[p]
        
        ## Updates 
        
        # Counter
        iteration_count <<- i*length(Approach)
        
        # Append the results data.frames
        res <- rbind(res, resIteration)
        t1 <- Sys.time()
        
        # Time
        time[i, 'Time'] <- as.numeric(difftime(t1, t0, units = 'secs'))
        i <- i+1
      }
      
      # Intermediate saves
      file_path <- paste0('Scenario2/res_', R, '_repetitions_', (n_noise[p] + n_var_inf), '_', sep, '_', repa, '.txt')
      write.table(res, file_path)
      res <- NULL
    }
  }
}
close(pb)
parallel::stopCluster(cl)

###########################################################################################################################################################################
###########################################################################################################################################################################

## Uploading the results for graphics and tables --------------------------------------------------------------------------------------------------------------------------
# To upload the results of scenario 2 we have to gather all the files corresponding to the same situation. Hence we make a distinction with the number of features first.

# All the files in scenario 2
files <- dir('Scenario2/')

# Get files related to the same number of features 
files_50 <- paste('Scenario2/', files[grep('res_50_repetitions_50', files)], sep = '')
files_100 <- paste('Scenario2/', files[grep('res_50_repetitions_100', files)], sep = '')
files_150 <- paste('Scenario2/', files[grep('res_50_repetitions_150', files)], sep = '')
files_200 <- paste('Scenario2/', files[grep('res_50_repetitions_200', files)], sep = '')

# Aggregating the files in one dataframe
res_50 <- data.frame(do.call('rbind', lapply(files_50, read.table)))
res_100 <- data.frame(do.call('rbind', lapply(files_100, read.table)))
res_150 <- data.frame(do.call('rbind', lapply(files_150, read.table)))
res_200 <- data.frame(do.call('rbind', lapply(files_200, read.table)))

## Aggregation ------------------------------------------------------------------------------------------------------------------------------------------------------------

# For 50 features
RES_50 <- data.frame(Approach = rep(Approach, each = length(separation)*length(repartition)), 
                     NbVar = res_50$NbVar[1],
                     Separation = rep(separation, each = length(repartition)), 
                     Repartition = rep(repartition, each = 1))

# For 100 features
RES_100 <- data.frame(Approach = rep(Approach, each = length(separation)*length(repartition)),
                      NbVar = res_100$NbVar[1],
                      Separation = rep(separation, each = length(repartition)),
                      Repartition = rep(repartition, each = 1))

# For 150 features
RES_150 <- data.frame(Approach = rep(Approach, each = length(separation)*length(repartition)),
                      NbVar = res_150$NbVar[1],
                      Separation = rep(separation, each = length(repartition)),
                      Repartition = rep(repartition, each = 1))

# For 200 features
RES_200 <- data.frame(Approach = rep(Approach, each = length(separation)*length(repartition)), 
                      NbVar = res_200$NbVar[1],
                      Separation = rep(separation, each = length(repartition)), 
                      Repartition = rep(repartition, each = 1))

## Filling the rows with mean results
for(i in 1:nrow(RES_50)){
  # For 50 features -------------------------------------
  lignes <- res_50[which(res_50$Approach == RES_50[i, 'Approach'] & 
                           res_50$NbVar == RES_50[i, 'NbVar'] & 
                           res_50$Separation == RES_50[i, 'Separation'] &
                           res_50$Repartition == RES_50[i, 'Repartition']), ]
  
  # Number of features
  RES_50[i, 'NbVarSelected'] <- mean(lignes[, 'NbVarSelected'])
  RES_50[i, 'NbVarInf'] <- mean(lignes[, 'NbVarInf'])
  
  # Selection indicators
  RES_50[i, 'Accuracy_selection'] <- mean(lignes[, 'Accuracy_selection'])*100
  RES_50[i, 'Specificity_selection'] <- mean(lignes[, 'Specificity_selection'])*100
  RES_50[i, 'Sensitivity_selection'] <- mean(lignes[, 'Sensitivity_selection'])*100
  RES_50[i, 'MCC_selection'] <- mean(lignes[, 'MCC_selection'])
  
  # Classification indicators (training data)
  RES_50[i, 'Accuracy_training'] <- mean(lignes[, 'Accuracy_training'])*100
  RES_50[i, 'Specificity_training'] <- mean(lignes[, 'Specificity_training'])*100
  RES_50[i, 'Sensitivity_training'] <- mean(lignes[, 'Sensitivity_training'])*100
  RES_50[i, 'MCC_training'] <- mean(lignes[, 'MCC_training'])
  
  # Classification indicators (test data)
  RES_50[i, 'Accuracy_test'] <- mean(lignes[, 'Accuracy_test'])*100
  RES_50[i, 'Specificity_test'] <- mean(lignes[, 'Specificity_test'])*100
  RES_50[i, 'Sensitivity_test'] <- mean(lignes[, 'Sensitivity_test'])*100
  RES_50[i, 'MCC_test'] <- mean(lignes[, 'MCC_test'])
  
  # For 100 features -------------------------------------
  lignes <- res_100[which(res_100$Approach == RES_100[i, 'Approach'] & 
                            res_100$NbVar == RES_100[i, 'NbVar'] & 
                            res_100$Separation == RES_100[i, 'Separation'] & 
                            res_100$Repartition == RES_100[i, 'Repartition']),]
  
  # Number of features
  RES_100[i, 'NbVarSelected'] <- mean(lignes[, 'NbVarSelected'])
  RES_100[i, 'NbVarInf'] <- mean(lignes[, 'NbVarInf'])
  
  # Selection indicators
  RES_100[i, 'Accuracy_selection'] <- mean(lignes[, 'Accuracy_selection'])*100
  RES_100[i, 'Specificity_selection'] <- mean(lignes[, 'Specificity_selection'])*100
  RES_100[i, 'Sensitivity_selection'] <- mean(lignes[, 'Sensitivity_selection'])*100
  RES_100[i, 'MCC_selection'] <- mean(lignes[, 'MCC_selection'])
  
  # Classification indicators (training data)
  RES_100[i, 'Accuracy_training'] <- mean(lignes[, 'Accuracy_training'])*100
  RES_100[i, 'Specificity_training'] <- mean(lignes[, 'Specificity_training'])*100
  RES_100[i, 'Sensitivity_training'] <- mean(lignes[, 'Sensitivity_training'])*100
  RES_100[i, 'MCC_training'] <- mean(lignes[, 'MCC_training'])
  
  # Classification indicators (test data)
  RES_100[i, 'Accuracy_test'] <- mean(lignes[, 'Accuracy_test'])*100
  RES_100[i, 'Specificity_test'] <- mean(lignes[, 'Specificity_test'])*100
  RES_100[i, 'Sensitivity_test'] <- mean(lignes[, 'Sensitivity_test'])*100
  RES_100[i, 'MCC_test'] <- mean(lignes[, 'MCC_test'])  
  
  # For 150 features -------------------------------------
  lignes <- res_150[which(res_150$Approach == RES_150[i, 'Approach'] & 
                            res_150$NbVar == RES_150[i, 'NbVar'] & 
                            res_150$Separation == RES_150[i, 'Separation'] & 
                            res_150$Repartition == RES_150[i, 'Repartition']),]
  
  # Number of features
  RES_150[i, 'NbVarSelected'] <- mean(lignes[, 'NbVarSelected'])
  RES_150[i, 'NbVarInf'] <- mean(lignes[, 'NbVarInf'])
  
  # Selection indicators
  RES_150[i, 'Accuracy_selection'] <- mean(lignes[, 'Accuracy_selection'])*100
  RES_150[i, 'Specificity_selection'] <- mean(lignes[, 'Specificity_selection'])*100
  RES_150[i, 'Sensitivity_selection'] <- mean(lignes[, 'Sensitivity_selection'])*100
  RES_150[i, 'MCC_selection'] <- mean(lignes[, 'MCC_selection'])
  
  # Classification indicators (training data)
  RES_150[i, 'Accuracy_training'] <- mean(lignes[, 'Accuracy_training'])*100
  RES_150[i, 'Specificity_training'] <- mean(lignes[, 'Specificity_training'])*100
  RES_150[i, 'Sensitivity_training'] <- mean(lignes[, 'Sensitivity_training'])*100
  RES_150[i, 'MCC_training'] <- mean(lignes[, 'MCC_training'])
  
  # Classification indicators (test data)
  RES_150[i, 'Accuracy_test'] <- mean(lignes[, 'Accuracy_test'])*100
  RES_150[i, 'Specificity_test'] <- mean(lignes[, 'Specificity_test'])*100
  RES_150[i, 'Sensitivity_test'] <- mean(lignes[, 'Sensitivity_test'])*100
  RES_150[i, 'MCC_test'] <- mean(lignes[, 'MCC_test'])  
  
  # For 200 features -------------------------------------
  lignes <- res_200[which(res_200$Approach == RES_200[i, 'Approach'] & 
                            res_200$NbVar == RES_200[i, 'NbVar'] & 
                            res_200$Separation == RES_200[i, 'Separation'] & 
                            res_200$Repartition == RES_200[i, 'Repartition']),]
  
  # Number of features
  RES_200[i, 'NbVarSelected'] <- mean(lignes[, 'NbVarSelected'])
  RES_200[i, 'NbVarInf'] <- mean(lignes[, 'NbVarInf'])
  
  # Selection indicators
  RES_200[i, 'Accuracy_selection'] <- mean(lignes[, 'Accuracy_selection'])*100
  RES_200[i, 'Specificity_selection'] <- mean(lignes[, 'Specificity_selection'])*100
  RES_200[i, 'Sensitivity_selection'] <- mean(lignes[, 'Sensitivity_selection'])*100
  RES_200[i, 'MCC_selection'] <- mean(lignes[, 'MCC_selection'])
  
  # Classification indicators (training data)
  RES_200[i, 'Accuracy_training'] <- mean(lignes[, 'Accuracy_training'])*100
  RES_200[i, 'Specificity_training'] <- mean(lignes[, 'Specificity_training'])*100
  RES_200[i, 'Sensitivity_training'] <- mean(lignes[, 'Sensitivity_training'])*100
  RES_200[i, 'MCC_training'] <- mean(lignes[, 'MCC_training'])
  
  # Classification indicators (test data)
  RES_200[i, 'Accuracy_test'] <- mean(lignes[, 'Accuracy_test'])*100
  RES_200[i, 'Specificity_test'] <- mean(lignes[, 'Specificity_test'])*100
  RES_200[i, 'Sensitivity_test'] <- mean(lignes[, 'Sensitivity_test'])*100
  RES_200[i, 'MCC_test'] <- mean(lignes[, 'MCC_test'])  
}

## Code for Table 3 --------------------------------------------------------------------------------------------------------------------------------------------------------

# Upper part of the table
print(xtable::xtable(cbind(RES_200[which(RES_200$Separation == 'forte' & RES_200$Repartition == 'equilibre'), 
              c('Approach', 'NbVarSelected', 'NbVarInf', 'Specificity_selection', 'Sensitivity_selection', 'MCC_test')],
      RES_200[which(RES_200$Separation == 'faible' & RES_200$Repartition == 'equilibre'), 
              c('NbVarSelected', 'NbVarInf', 'Specificity_selection', 'Sensitivity_selection', 'MCC_test')]), style = 'latex', digits = c(0, 0, 2, 2, 2, 2, 3, 2, 2, 2, 2, 3)))

# Lower part of the table
print(xtable::xtable(cbind(RES_200[which(RES_200$Separation == 'forte' & RES_200$Repartition == 'desequilibre'),
              c('Approach', 'NbVarSelected', 'NbVarInf', 'Specificity_selection', 'Sensitivity_selection', 'MCC_test')],
      RES_200[which(RES_200$Separation == 'faible' & RES_200$Repartition == 'desequilibre'),
              c('NbVarSelected', 'NbVarInf', 'Specificity_selection', 'Sensitivity_selection', 'MCC_test')]), style = 'latex', digits = c(0, 0, 2, 2, 2, 2, 3, 2, 2, 2, 2, 3)))





## Code for Figure 3 -----------------------------------------------------------------------------------------------------------------------------------------------------
# To plot the heatmap of the selected features we build a dataframe on a long format with how much time each feature is selected

# Dataframe for each situation for 200 features
VARSELECTED_LONG_faible_desequilibre <- dataframeSelectedFeaturesScenario2(res_200[which(res_200$Separation == 'faible' & res_200$Repartition == 'desequilibre'),])
VARSELECTED_LONG_faible_equilibre <- dataframeSelectedFeaturesScenario2(res_200[which(res_200$Separation == 'faible' & res_200$Repartition == 'equilibre'),])
VARSELECTED_LONG_forte_desequilibre <- dataframeSelectedFeaturesScenario2(res_200[which(res_200$Separation == 'forte' & res_200$Repartition == 'desequilibre'),])
VARSELECTED_LONG_forte_equilibre <- dataframeSelectedFeaturesScenario2(res_200[which(res_200$Separation == 'forte' & res_200$Repartition == 'equilibre'),])

# Modification of the labels 
VARSELECTED_LONG_faible_desequilibre[, c('Balance', 'Separation')] <- rep(c('Unbalanced', 'Weak'), each = nrow(VARSELECTED_LONG_faible_desequilibre))
VARSELECTED_LONG_faible_equilibre[, c('Balance', 'Separation')] <- rep(c('Balanced', 'Weak'), each = nrow(VARSELECTED_LONG_faible_equilibre))
VARSELECTED_LONG_forte_desequilibre[, c('Balance', 'Separation')] <- rep(c('Unbalanced', 'Strong'), each = nrow(VARSELECTED_LONG_forte_desequilibre))
VARSELECTED_LONG_forte_equilibre[, c('Balance', 'Separation')] <- rep(c('Balanced', 'Strong'), each = nrow(VARSELECTED_LONG_forte_equilibre))

# Binding to a single dataframe for plot
VARSELECTED_LONG <- rbind(VARSELECTED_LONG_faible_desequilibre, VARSELECTED_LONG_faible_equilibre, VARSELECTED_LONG_forte_desequilibre, VARSELECTED_LONG_forte_equilibre)

# Plot 
plot_figure_03 <- ggplot(data = VARSELECTED_LONG, aes(x = Feature, y = variable))+
  geom_tile(aes(fill = value))+
  scale_x_discrete(breaks = unique(VARSELECTED_LONG$Feature), limits = unique(VARSELECTED_LONG$Feature))+
  scale_y_discrete(limits = rev(Approach))+
  scale_fill_gradientn(colours = rev(viridis::inferno(10)), breaks = seq(0, R, R/5))+
  labs(fill = '', x = '', y =  '')+
  facet_grid(Balance~Separation)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 3, hjust = 1),
        axis.text.y = element_text(size = 20),
        legend.key.width = unit(0.75, 'cm'),
        legend.key.height = unit(2, 'cm'),
        legend.margin = margin(0.1, 0.1, 0.1, 0.1),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 30, margin = margin(0.025, 0, 0.025, 0, 'cm')),
        strip.text.y = element_text(size = 30, margin = margin(0, 0.025, 0, 0.025, 'cm')))

# Save the plots
setEPS()
postscript('Figures/Figure_03.eps', width = 18, height = 10, horizontal = FALSE)
plot_figure_03
dev.off()

## Code for Figure 4 --------------------------------------------------------------------------------------------------------------------------------------------------------

# Changing the labels
res_200$Separation <- factor(res_200$Separation, levels = c('forte', 'faible'), labels = c('Strong', 'Weak'))
res_200$Repartition <- factor(res_200$Repartition, levels = c('equilibre', 'desequilibre'), labels = c('Balanced', 'Unbalanced'))

# Plot
plot_figure_04 <- ggplot(data = res_200, aes(x = Approach, y = MCC_test))+
  geom_boxplot(size = 0.3, outlier.size = 0.8, outlier.shape = 19)+
  coord_flip()+
  scale_x_discrete(limits = rev(Approach), labels = rev(Approach))+
  scale_y_continuous(limits = c(-0.5, 1))+
  grids(axis = 'y', color = 'grey', linetype = 'dashed', size = 0.2)+
  grids(axis = 'x', color = 'grey', linetype = 'solid', size = 0.1)+
  labs(x = '', y = '')+
  facet_grid(Repartition~Separation)+
  theme(axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 20),
        panel.background = element_rect(fill = 'white', colour = 'black', linewidth = 0.2),
        strip.text.x = element_text(size = 25, margin = margin(0.05, 0, 0.05, 0, 'cm')), 
        strip.text.y = element_text(size = 25, margin = margin(0, 0.05, 0, 0.05, 'cm')))

setEPS()
postscript('Figures/Figure_04.eps', width = 18, height = 10, horizontal = FALSE)
plot_figure_04
dev.off()
