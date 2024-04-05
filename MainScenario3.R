####################################################################################################################################
#                                                           SCENARIO 3                                                             #  
####################################################################################################################################

## Environment ---------------------------------------------------------------------------------------------------------------------
if(!require(devtools)){
  # Package to install package
  install.packages('devtools')
}

if(!require(stringr)){
  # Package to manipulate string
  devtools::install_version('stringr', version = '1.5.0')
  library(stringr)
}

if(!require(parallel)){
  # Package to execute parallel work
  devtools::install_version('parallel', version = '4.3.1')
  library(parallel)
}

if(!require(doSNOW)){
  # Package to execute parallel work
  devtools::install_version('doSNOW', version = '1.0.20')
  library(doSNOW)
}

if(!require(ggplot2)){
  # Package use for plots
  devtools::install_version('ggplot2', version = '3.4.3')
  library(ggplot2)
}

if(!require(ggpubr)){
  # Package used for grids in ggplots
  devtools::install_version('ggpubr', version = '0.6.0')
  library(ggpubr)
}

## Sources ------------------------------------------------------------------------------------------------------------------------
source('functionGammaMetric.R')
source('functionGeneration.R')
source('functions.R')

## Parameters of the generation ---------------------------------------------------------------------------------------------------
N <- 2000                                                                   # Number of observations
s_g <- 10                                                                   # Size of each groups
n_g <- 10                                                                   # Number of groups
non_zero_coeff <- 1                                                         # Informative features are the first feature of each group 
non_zero_group <- c(1, 2, 3, 4, 5, 10)                                      # Groups with non-zero coefficients
value <- 1.5                                                                # Beta coefficients value (when different from 0)
n_ind <- 1                                                                  # Number of independent group
simulation <- 1:6                                                           # Cases for this scenario 
R <- 50                                                                     # Number of repetitions
set.seed(123)                                                               # Seed
beta_vec <- produceBeta(s_g, n_g, non_zero_coeff, non_zero_group, value)    # Beta vector
beta <- beta_vec[-1]
beta0 <- beta_vec[1]


## Parameters of the parallel computations ---------------------------------------------------------------------------------------
# Feature selection methods 
Approach <- c('BASELINE', 'BEST', 'CFS', 'CHI2', 'CONS', 'IG', 'IGR', 'ONER', 'RELIEF', 'RFI', 'SU', 'SVM-RFE', 
              'GAMMA_BACK', 'GAMMA_BF', 'GAMMA_FORW', 'GAMMA_HC')

# Number of clusters
n_clusters <- parallel::detectCores() - 1                             # Number of clusters available
cl <- parallel::makeCluster(n_clusters)                               # Initializing n_clusters
doSNOW::registerDoSNOW(cl)                                            # Parallel work

# Time vector
time <- expand.grid(Repetition = 1:R, 
                    Simulation = simulation, 
                    Time = 0)

# Initialisations
iteration_count <<- 0                                                 # Counter
i <- 1                                                                # Counter for progress bar
res <- NULL                                                           # Matrix of results

# Progress bar
total_iteration <- R*length(simulation)*length(Approach)              # Total number of iterations
pb <- txtProgressBar(min = 0, max = total_iteration, style = 3)       # A progress bar
progress <- function(n) setTxtProgressBar(pb, n + iteration_count)    # Function to update the progress bar
opts <- list(progress = progress)                                     # List of options for the progress bar

## BEGINING OF THE SIMULATION -----------------------------------------------------------------------------------------------------

for(simu in simulation){
  # We defined a set of parameters for generation of data with regard to the level and type of dependence we considered
  
  if(simu == 1){
    # First case with a constant correlation level alpha_max = 0.9
    alpha_max <- 0.9
    constant_cor <- TRUE
    
  }else if(simu == 2){
    # Second case with a constant correlation level alpha_max = 0.6
    alpha_max <- 0.6
    constant_cor <- TRUE
    
  }else if(simu == 3){
    # Third case with a constant correlation level alpha_max = 0.3
    alpha_max <- 0.3
    constant_cor <- TRUE
    
  }else if(simu == 4){
    # Fourth case with a non-constant correlation level:
    # Maximum level alpha_max = 0.9 and minimum level c = 0.35
    alpha_max <- 0.9
    c <- 0.35
    constant_cor <- FALSE
    
  }else if(simu == 5){
    # Fith case with a non-constant correlation level:
    alpha_max <- 0.6
    c <- 0.25
    constant_cor <- FALSE
  }else if(simu == 6){
    # Sixth case with a non-constant correlation level:
    alpha_max <- 0.3
    c <- 0.1
    constant_cor <- FALSE
  }
  
  ## Beginning of the repetitions
  for(r in 1:R){
    # Time counter
    t0 <- Sys.time()
    
    # Generation of the data (train data set)
    dat_train <- genData(n = N, s_g = s_g, n_g = n_g, beta0 = beta0, beta = beta, alpha_max = alpha_max, c = c, n_ind = n_ind, constant_cor = constant_cor)
    X_train <- dat_train[, -1]
    Y_train <- dat_train[, 1]
    
    # Generation of the data (test data set)
    dat_test <- genData(n = N, s_g = s_g, n_g = n_g, beta0 = beta0, beta = beta, alpha_max = alpha_max, c = c, n_ind = n_ind, constant_cor = constant_cor)
    X_test <- dat_test[, -1]
    Y_test <- dat_test[, 1]
    
    # Perform feature selection, modelisation and prediction in parallel for each feature selection method
    resIteration <- foreach(fs = Approach, .combine = 'rbind', .export = ls()[!ls() %in% c('performSimulationIteration', 'X_train', 'X_test', 'Y_train', 'Y_test', 'beta', 'c')], .options.snow = opts, .verbose = FALSE) %dopar%{
      localRes <- performSimulationIteration(X_train, Y_train, X_test, Y_test, beta, fs)
      row.names(localRes) <- NULL
      return(localRes)
    }
    
    # Filling the common columns
    resIteration$Simulation <- simu
    
    # Updates
    t1 <- Sys.time()
    iteration_count <- i*length(Approach)
    res <- rbind(res, resIteration)
    time[i, 'Time'] <- as.numeric(difftime(t1, t0, units = 'secs'))
    i <- i+1
  }
  
  # Intermediate save
  file_path <- paste0('Scenario3/res_simulation_', simu, '_', R, '_repetitions_', Sys.Date(), '.txt')
  write.table(res, file_path)
  res <- NULL
}
close(pb)
parallel::stopCluster(cl)

###########################################################################################################################################################################
###########################################################################################################################################################################

## Uploading the results for graphics and tables --------------------------------------------------------------------------------------------------------------------------

# All the results files
files <- dir('Scenario3/')

# Aggregations of the results files
patterns <- paste0('res_simulation_', 1:6)
res <- data.frame(do.call('rbind', lapply(
                 paste0('Scenario3/', files[grep(paste(patterns, collapse = '|'), files)]), read.table
                 )
               ))

## Aggregation ------------------------------------------------------------------------------------------------------------------------------------------------------------

RES <- data.frame(expand.grid(
  Approach = Approach, 
  Simulation = simulation
))

for(i in 1:nrow(RES)){
  # Line
  lignes <- res[which(res$Approach == RES[i, 'Approach'] & res$Simulation == RES[i, 'Simulation']),]
  
  # Number of features
  RES[i, 'NbVarSelected'] <- mean(lignes[,'NbVarSelected'])
  RES[i, 'NbVarInf'] <- mean(lignes[, 'NbVarInf'])
  
  # Selection indicators
  RES[i, 'Accuracy_selection'] <- mean(lignes[,'Accuracy_selection'])*100
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

#### Table ----------------------------------------------------------------------------------------------------------------------------------------------------------------

## Code for Table 4 -------------------------------------------------------------------------------------------------------------------------------------------------------

# Upper part of the table
cbind(RES[which(RES$Simulation == 1), c('Approach', 'NbVarSelected', 'NbVarInf', 'Specificity_selection', 'Sensitivity_selection', 'Accuracy_test')],
      RES[which(RES$Simulation == 4), c('NbVarSelected', 'NbVarInf', 'Specificity_selection', 'Sensitivity_selection', 'Accuracy_test')])

# Middle part of the table
cbind(RES[which(RES$Simulation == 2), c('Approach', 'NbVarSelected', 'NbVarInf', 'Specificity_selection', 'Sensitivity_selection', 'Accuracy_test')],
      RES[which(RES$Simulation == 5), c('NbVarSelected', 'NbVarInf', 'Specificity_selection', 'Sensitivity_selection', 'Accuracy_test')])

# Lower part of the table
cbind(RES[which(RES$Simulation == 3), c('Approach', 'NbVarSelected', 'NbVarInf', 'Specificity_selection', 'Sensitivity_selection', 'Accuracy_test')],
      RES[which(RES$Simulation == 4), c('NbVarSelected', 'NbVarInf', 'Specificity_selection', 'Sensitivity_selection', 'Accuracy_test')])

#### Figures --------------------------------------------------------------------------------------------------------------------------------------------------------------

## Code for Figure 5 ------------------------------------------------------------------------------------------------------------------------------------------------------
# Dataframe for each simulation situation
VARSELECTED_LONG_constant_high <- dataframeSelectedFeaturesScenario3(res[res$Simulation == 1,])
VARSELECTED_LONG_constant_med <- dataframeSelectedFeaturesScenario3(res[res$Simulation == 2,])
VARSELECTED_LONG_constant_low <- dataframeSelectedFeaturesScenario3(res[res$Simulation == 3,])

VARSELECTED_LONG_non_constant_high <- dataframeSelectedFeaturesScenario3(res[res$Simulation == 4,])
VARSELECTED_LONG_non_constant_med <- dataframeSelectedFeaturesScenario3(res[res$Simulation == 5,])
VARSELECTED_LONG_non_constant_low <- dataframeSelectedFeaturesScenario3(res[res$Simulation == 6,])

# Modification of the labels 
VARSELECTED_LONG_constant_high[, c('Constant', 'Correlation')] <- rep(c(TRUE, '0.9'), each = nrow(VARSELECTED_LONG_constant_high))
VARSELECTED_LONG_constant_med[, c('Constant', 'Correlation')] <- rep(c(TRUE, '0.6'), each = nrow(VARSELECTED_LONG_constant_med))
VARSELECTED_LONG_constant_low[, c('Constant', 'Correlation')] <- rep(c(TRUE, '0.3'), each = nrow(VARSELECTED_LONG_constant_low))
VARSELECTED_LONG_non_constant_high[, c('Constant', 'Correlation')] <- rep(c(FALSE, '0.9'), each = nrow(VARSELECTED_LONG_non_constant_high))
VARSELECTED_LONG_non_constant_med[, c('Constant', 'Correlation')] <- rep(c(FALSE, '0.6'), each = nrow(VARSELECTED_LONG_non_constant_med))
VARSELECTED_LONG_non_constant_low[, c('Constant', 'Correlation')] <- rep(c(FALSE, '0.3'), each = nrow(VARSELECTED_LONG_non_constant_low))

# Merging into one dataframe
VARSELECTED_LONG <- rbind(
  VARSELECTED_LONG_constant_high, VARSELECTED_LONG_constant_med, VARSELECTED_LONG_constant_low, 
  VARSELECTED_LONG_non_constant_high, VARSELECTED_LONG_non_constant_med, VARSELECTED_LONG_non_constant_low
)

# Some formatting for plot aesthetics
VARSELECTED_LONG$Constant <- factor(VARSELECTED_LONG$Constant, 
                                    levels = c(TRUE, FALSE), 
                                    labels = c('CONSTANT', 'NON-CONSTANT'))

VARSELECTED_LONG$Correlation <- factor(VARSELECTED_LONG$Correlation, 
                                       levels = c('0.9', '0.6', '0.3'),
                                       labels = c(
                                         expression(paste(alpha[max], ' = 0.9')),
                                         expression(paste(alpha[max], ' = 0.6')),
                                         expression(paste(alpha[max], ' = 0.3'))
                                       ))
# X-axis labels to display non-informative features as N with the index of the feature
index_of_informative_features <- s_g*(non_zero_group - 1) + non_zero_coeff                  # It only works if non_zero_coeff is a single element vector
x_labels <- paste0(1:(s_g*n_g))
x_labels[index_of_informative_features] <- paste0('x', x_labels[index_of_informative_features])
x_labels[-index_of_informative_features] <- paste0('N', x_labels[-index_of_informative_features])

# Figure
plot_selection <- ggplot(data = VARSELECTED_LONG, aes(x = Feature, y = variable))+
  geom_tile(aes(fill = value), color = 'grey', lwd = 0.05, linetype = 1)+
  scale_x_discrete(breaks = unique(VARSELECTED_LONG$Feature), limits = unique(VARSELECTED_LONG$Feature), labels = x_labels)+
  scale_y_discrete(limits = rev(Approach))+
  scale_fill_gradientn(colours = rev(viridis::inferno(10)), breaks = seq(0, R, R/5))+
  labs(fill = '', x = '', y = '')+
  facet_grid(Correlation~Constant, labeller = label_parsed)+
  theme(axis.text.x = element_text(angle = 90, size = 2, vjust = 0.5, hjust = 0),
        axis.text.y = element_text(size = 4),
        legend.key.width = unit(0.25, 'cm'),
        legend.key.size = unit(0.2, 'cm'),
        legend.key.height = unit(0.5, 'cm'),
        legend.margin = margin(0.01, 0.01, 0.01, 0.01),
        legend.text = element_text(size = 5),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.length = unit(0.025, 'cm'),
        strip.text.x = element_text(size = 5, margin = margin(0.025, 0, 0.025, 0, 'cm')),
        strip.text.y = element_text(size = 5, margin = margin(0, 0.025, 0, 0.025, 'cm')))

# Saving the plot
tiff(filename = 'temp_res/Scenario3/Ngo_Figure5.tiff', width = 1750*2, height = 936*2, res = 600, units = 'px', pointsize = 12)
plot_selection
dev.off()


## Code for Figure 6 ------------------------------------------------------------------------------------------------------------------------------------------------------

# Changing the labels
res[which(res$Simulation == 1), c('Constant', 'Correlation')] <- rep(c(TRUE, '0.9'), each = nrow(res[which(res$Simulation == 1),]))
res[which(res$Simulation == 2), c('Constant', 'Correlation')] <- rep(c(TRUE, '0.6'), each = nrow(res[which(res$Simulation == 2),]))
res[which(res$Simulation == 3), c('Constant', 'Correlation')] <- rep(c(TRUE, '0.3'), each = nrow(res[which(res$Simulation == 3),]))

res[which(res$Simulation == 4), c('Constant', 'Correlation')] <- rep(c(FALSE, '0.9'), each = nrow(res[which(res$Simulation == 4),]))
res[which(res$Simulation == 5), c('Constant', 'Correlation')] <- rep(c(FALSE, '0.6'), each = nrow(res[which(res$Simulation == 5),]))
res[which(res$Simulation == 6), c('Constant', 'Correlation')] <- rep(c(FALSE, '0.3'), each = nrow(res[which(res$Simulation == 6),]))

# Constant and Correlation as factors•
res$Constant <- factor(res$Constant, levels = c(TRUE, FALSE), labels = c('CONSTANT', 'NON-CONSTANT'))
res$Correlation <- factor(res$Correlation, 
                          levels = c('0.9', '0.6', '0.3'),
                          labels = c(
                            expression(paste(alpha[max], ' = 0.9')),
                            expression(paste(alpha[max], ' = 0.6')),
                            expression(paste(alpha[max], ' = 0.3'))
                          ))

# Plot
plot_accuracy <- ggplot(data = res, aes(x = Approach, y = Accuracy_test))+
  geom_boxplot(size = 0.2, outlier.size = 0.05, outlier.shape = 19)+
  coord_flip()+
  scale_x_discrete(limits = rev(Approach), labels = rev(Approach))+
  scale_y_continuous(limits = c(0.45, 1))+
  grids(axis = 'y', color = 'grey', linetype = 'dashed', size = 0.2)+
  grids(axis = 'x', color = 'grey', linetype = 'solid', size = 0.1)+
  labs(x = '', y = '')+
  facet_grid(Correlation~Constant, labeller = label_parsed)+
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 4),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.length = unit(0.025, 'cm'),
        panel.background = element_rect('white', colour = 'black', linewidth = 0.2),
        strip.text.x = element_text(size = 7, margin = margin(0.05, 0, 0.05, 0, 'cm')),
        strip.text.y = element_text(size = 7, margin = margin(0, 0.05, 0, 0.05, 'cm')))

# Saving plots
tiff(filename = 'temp_res/Scenario3/Ngo_Figure6.tiff', width = 1750*2, height = 936*2, res = 600, units = 'px', pointsize = 12)
plot_accuracy
dev.off()



