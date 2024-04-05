###########################################################################################################
#                                                 SCENARIO 1                                              #
###########################################################################################################

## Environment --------------------------------------------------------------------------------------------
if(!require(devtools)){
  # Package used to install packages
  install.package(devtools)
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

if(!require(ggpubr)){
  # Package used for grids in plots
  devtools::install_version('ggpubr', version = '0.6.0')
  library(ggpubr)
}
#### Sources --------------------------------------------------------------------------------------------
source('functionGammaMetric.R')
source('functionGeneration.R')
source('functions.R')

#### Functions ------------------------------------------------------------------------------------------

# List of functions with regard to the method
# approcheFunction <- list(
#   'BASELINE' <- function(f, data) colnames(data[, -1]),
#   'BEST' <- function(f, data) colnames(data[, -1])[1:n_var_inf], 
#   'CFS' <- function(f, data) FSelector::cfs(f, data), 
#   'CHI2' <- function(f, data) FSelector::cutoff.biggest.diff(FSelector::chi.squared(f, data)),
#   'CONS' <- function(f, data) FSelector::consistency(f, data),
#   'IG' <- function(f, data) FSelector::cutoff.biggest.diff(FSelector::information.gain(f, data)),
#   'IGR' <- function(f, data) FSelector::cutoff.biggest.diff(FSelector::gain.ratio(f, data)),
#   'ONER' <- function(f, data) FSelector::cutoff.biggest.diff(FSelector::oneR(f, data)),
#   'RELIEF' <- function(f, data) FSelector::cutoff.biggest.diff(FSelector::relief(f, data)),
#   'RFI' <- function(f, data) FSelector::cutoff.biggest.diff(FSelector::random.forest.importance(f, data)),
#   'SU' <- function(f, data) FSelector::cutoff.biggest.diff(FSelector::symmetrical.uncertainty(f, data)),
#   'SVM-RFE' <- function(f, data){
#     if(!require('mlr3verse')){
#       # Package used to operate SVM-RFE
#       devtools::install_version('mlr3verse', version = '0.2.8')
#       library(mlr3verse)
#     }
#     
#     # We use a Recursive Feature Elimination with a SVM classifier
#     
#     # Just a checking
#     data$y <- factor(data$y)
#     
#     # Definition 
#     optimizer <- mlr3verse::fs('rfe', n_features = 1, feature_number = 1)
#     
#     # Definition of the learner
#     learner <- mlr3::lrn('classif.svm', type = 'C-classification', kernel = 'linear', predict_type = 'prob')
#     
#     # Definition of the task
#     task <- mlr3::TaskClassif$new(id = 'simulated_data', backend = data, target = 'y', positive = '1')
#     
#     # Add the task to the dictionary
#     mlr_tasks$add('simulated_data', task)
#     
#     # Instance
#     instance <- mlr3verse::fsi(task = tsk('simulated_data'), learner = learner, resampling = rsmp('cv', folds = 5), 
#                                measures = msr('classif.mcc'), terminator = trm('none'), callback = clbk('mlr3fselect.svm_rfe'))
#     
#     # Optimizer
#     optimizer$optimize(inst = instance)
#     var_select <- instance$result_feature_set
#     
#     return(var_select)
#   }, 
#   'GAMMA_BACK' <- function(f, data) FSelector::backward.search(colnames(data[, -1]), eval.fun = function(subset){
#     gammaMetric(X_train[, subset], class = Y_train, covEstim = 'empiric', plot = FALSE)
#   }),
#   'GAMMA_BF' <- function(f, data) FSelector::best.first.search(colnames(data[, -1]), eval.fun = function(subset){
#     gammaMetric(X_train[, subset], class = Y_train, covEstim = 'empiric', plot = FALSE)
#   }),
#   'GAMMA_FORW' <- function(f, data) FSelector::forward.search(colnames(data[, -1]), eval.fun = function(subset){
#     gammaMetric(X_train[, subset], class = Y_train, covEstim = 'empiric', plot = FALSE)
#   }),
#   'GAMMA_HC' <- function(f, data) FSelector::hill.climbing.search(colnames(data[, -1]), eval.fun = function(subset){
#     gammaMetric(X_train[, subset], class = Y_train, covEstim = 'empiric', plot = FALSE)
#   })
# )

# # Iteration function
# performSimulationIteration <- function(X_train, Y_train, X_test, Y_test, fs){
#   # Feature selection step
#   var_select <- do.call(fs, list(f = FSelector::as.simple.formula(colnames(X_train), 'y'), data = data.frame(y = factor(Y_train), X_train)))
#   
#   # Type of features selected 
#   informative_features <- grep(pattern = 'x', x = colnames(X_train), value = TRUE)
#   selected <- (colnames(X_train) %in% var_select)*1
#   informative <- (colnames(X_train) %in% informative_features)*1
#   
#   # Correctly and uncorrectly selected features
#   cm_selection <- caret::confusionMatrix(data = factor(selected, levels = c('0', '1')), reference = factor(informative), positive = '1')
#   acc_sel <- cm_selection$overall['Accuracy']
#   spe_sel <- cm_selection$byClass['Specificity']
#   sen_sel <- cm_selection$byClass['Sensitivity']
#   mcc_sel <- mltools::mcc(preds = factor(selected, levels = c('0', '1')), actuals = factor(informative))
#   
#   # Construction of the logistic model with the selected features
#   ctr <- caret::trainControl(method = 'none', allowParallel = FALSE)
#   mod.fit <- caret::train(FSelector::as.simple.formula(var_select, 'y'), 
#                           data = data.frame(y = factor(Y_train), X_train), 
#                           trControl = ctr, 
#                           method = 'glm', 
#                           family = 'binomial')
#   
#   # Get the estimate beta coefficients
#   resume <- summary(mod.fit)
#   coef_estimate <- resume$coef[, 'Estimate']
#   std_estimate <- resume$coef[, 'Std. Error']
#   
#   # Create vector with NA for non selected features and beta estimates (and std of the estimate)
#   estimates <- data.frame(Beta = rep(NA, length(colnames(X_train))),
#                           Std = rep(NA, length(colnames(X_train))))
#   row.names(estimates) <- colnames(X_train)
#   
#   estimates[names(coef_estimate)[-1], c('Beta', 'Std')] <- c(coef_estimate[-1], std_estimate[-1])
#   
#   # Fitted values of the training sample
#   fitted <- predict(mod.fit, data.frame(y = Y_train, X_train))
#   
#   # Classification performances (training sample)
#   cm <- caret::confusionMatrix(fitted, factor(Y_train), positive = '1')
#   acc_train <- cm$overall['Accuracy']
#   spe_train <- cm$byClass['Specificity']
#   sen_train <- cm$byClass['Sensitivity']
#   mcc_train <- mltools::mcc(fitted, factor(Y_train))
#   
#   # Prediction and classification performances (test sample)
#   preds <- predict(mod.fit, data.frame(y = Y_test, X_test))
#   cm <- caret::confusionMatrix(preds, factor(Y_test), positive = '1')
#   acc_test <- cm$overall['Accuracy']
#   spe_test <- cm$byClass['Specificity']
#   sen_test <- cm$byClass['Sensitivity']
#   mcc_test <- mltools::mcc(preds, factor(Y_test))
#   
#   # Aggregation of the results in a dataframe
#   localRes <- data.frame(Approach = fs, NbVarSelected = length(var_select), NbVarInf = length(grep('x', var_select)), 
#                          Accuracy_selection = acc_sel, Specificity_selection = spe_sel, Sensitivity_selection = sen_sel, MCC_selection = mcc_sel, 
#                          Accuracy_training = acc_train, Specificity_training = spe_train, Sensitivity_training = sen_train, MCC_training = mcc_train, 
#                          Accuracy_test = acc_test, Specificity_test = spe_test, Sensitivity_test = sen_test, MCC_test = mcc_test, 
#                          VarSelected = paste(var_select, collapse = ';'), coef_estimate = paste0(estimates$Beta, collapse = ';'), std_estimate = paste0(estimates$Std, collapse = ';'))
#   
#   return(localRes)
# }

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
  resIteration <- foreach(fs = Approach, .combine = 'rbind', .export = ls()[!ls() %in% c('performSimulationIteration', 'X_train', 'Y_train', 'X_test', 'Y_test')], .options.snow = opts, .verbose = FALSE) %dopar% {
    localRes <- performSimulationIteration_draft(X_train, Y_train, X_test, Y_test, beta = beta, fs)
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
file_path <- paste0('Scenario1/res_', R, '_iterations_', Sys.Date(), '.txt')
write.table(res, file_path)

############################################################################################################################################################################
############################################################################################################################################################################

#### Uploading the results for graphics and tables -------------------------------------------------------
res <- read.table('Scenario1/res_50_iterations_2024-02-20.txt')

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
Selection_plot <- ggplot(data = VARSELECTED_LONG, aes(x = Feature, y = variable))+
  geom_tile(aes(fill = value), color = 'white', lwd = 0.3, linetype = 1)+
  scale_x_discrete(breaks = c(paste0('x', 1:n_var_inf), paste0('Noise', 1:n_noise)), limits = c(paste0('x', 1:n_var_inf), paste0('Noise', 1:n_noise)), labels = c(paste0('x', 1:n_var_inf), paste0('N', 1:n_noise)))+
  scale_y_discrete(limits = rev(Approach))+
  scale_fill_gradientn(colours = rev(viridis::inferno(10)), breaks = seq(0, R, R/5))+
  labs(fill = '', x = '', y = '')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, size = 10),
        axis.title = element_text(size= 7), 
        title = element_text(size = 7),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.length = unit(0.025, 'cm'))

tiff(filename = 'Figures/Ngo_Figure1.tiff', width = 2*1750, height = 2*936, units = 'px', pointsize = 10, res = 600)
Selection_plot
dev.off()

## Code for Figure 2 -------------------------------------------------------------------------------------
### Plot for classification performances
plot_acc <- ggplot(data = res, aes(x = Approach, y = Accuracy_test))+
  geom_boxplot(size = 0.3, outlier.size = 0.02, outlier.shape = 19) +
  coord_flip()+
  scale_x_discrete(limits = rev(Approach), labels = rev(Approach))+
  scale_y_continuous(limits = c(0.5, 1))+
  grids(axis = 'y', color = 'grey', linetype = 'dashed', size = 0.2)+
  grids(axis = 'x', color = 'grey', linetype = 'solid', size = 0.2)+
  labs(y = '', x = '')+
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 7),
        title = element_text(size = 7),
        panel.background = element_rect(fill = 'white', colour = 'black', linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.length = unit(0.025, 'cm'))

tiff(filename = 'Figures/Ngo_Figure2.tiff', width = 2*1750, height = 2*936, units = 'px', pointsize = 10, res = 600)
plot_acc
dev.off()
  

####
# BONUS
####

res_new <- res
res <- read.table('Scenario1/res_50_iterations_2024-02-20.txt')

diff <- res_new[, c('NbVarSelected', 'NbVarInf', 'Accuracy_selection', 'Accuracy_test')] - res[, c('NbVarSelected', 'NbVarInf', 'Accuracy_selection', 'Accuracy_validation')]
diff$Approach <- res$Approach
diff <- diff[, c(5, 1:4)]
diff$Set <- rep(1:50, each = 16)
diff[which(diff$NbVarSelected != 0), ]
