#########################################################################################################################################
###                                                        SET OF FUNCTIONS FOR THE SIMULATION                                        ###
#########################################################################################################################################

## Packages -----------------------------------------------------------------------------------------------------------------------------
if(!require(FSelector)){
  # Package used for feature selection
  devtools::install_version('FSelector', version = 'FSelctor', version = 0.34)
  library(FSelector)
}

if(!require(caret)){
  # Package used to calibrate models (and compute some indicators)
  devtools::install_version('caret', version = '6.0.94')
  library(caret)
}

if(!require(reshape2)){
  # Package use to format data table
  devtools::install_version('reshape2', version = '1.4.4')
}

## Sources ------------------------------------------------------------------------------------------------------------------------------
source('functionGammaMetric.R')
source('functionGeneration.R')

## Feature selection methods ------------------------------------------------------------------------------------------------------------
# To facilitate parallel work, we put all the feature selection methods function into a list and R we go through that list to execute a 
# method and call the corresponding function. All the functions in that list have the same inputs and outputs. 
# Inputs:
#      f: The formula object we want to do selection for.
#   data: The dataset with the first column being the target feature y
# Outputs:

approachFunctions <- list(
  'BASELINE' <- function(f, data) colnames(data[, -1]),
  'BEST' <- function(f, data) colnames(data[, which(beta != 0)+1]),
  'CFS' <- function(f, data) FSelector::cfs(f, data),
  'CHI2' <- function(f, data) FSelector::cutoff.biggest.diff(FSelector::chi.squared(f, data)),
  'CONS' <- function(f, data) FSelector::consistency(f, data),
  'IG' <- function(f, data) FSelector::cutoff.biggest.diff(FSelector::information.gain(f, data)),
  'IGR' <- function(f, data) FSelector::cutoff.biggest.diff(FSelector::gain.ratio(f, data)),
  'ONER' <- function(f, data) FSelector::cutoff.biggest.diff(FSelector::oneR(f, data)),
  'RFI' <- function(f, data) FSelector::cutoff.biggest.diff(FSelector::random.forest.importance(f, data)),
  'RELIEF' <- function(f, data) FSelector::cutoff.biggest.diff(FSelector::relief(f, data)),
  'SU' <- function(f, data) FSelector::cutoff.biggest.diff(FSelector::symmetrical.uncertainty(f, data)),
  'SVM-RFE' <- function(f, data){
    if(!require(mlr3verse)){
      devtools::install_version('mlr3verse')
      library(mlr3verse)
    }
    # To ensure that y is a factor
    data$y <- factor(data$y)
    
    # Definition of the feature selection
    optimizer <- mlr3verse::fs('rfe', n_features = 1, feature_number = 1)
    
    # Definition of the learner
    learner <- mlr3::lrn("classif.svm", type = 'C-classification', kernel = 'linear', predict_type = 'prob')
    
    # Definition of the task
    task <- mlr3::TaskClassif$new(id = 'simulated_data', 
                                  backend = data, 
                                  target = 'y', 
                                  positive = '1')
    
    # Add the task to the dictionary
    mlr_tasks$add('simulated_data', task)
    
    # Instance
    instance <- mlr3verse::fsi(task = tsk('simulated_data'),
                               learner = learner, 
                               resampling = rsmp('cv', folds = 5),
                               measures = msr('classif.mcc'),
                               terminator = trm('none'),
                               callback = clbk('mlr3fselect.svm_rfe'))
    
    # Optimizer
    optimizer$optimize(inst = instance)
    var_select <- instance$result_feature_set
    
    return(var_select)
  },
  'GAMMA_BF' <- function(f, data){
    # Which estimator of the variance-covariance matrix to use
    estim <- ifelse(min(table(data[, 1])) < ncol(data[, -1]), 'shrink', 'empiric')
    
    # Feature selection with the best-first search direction
    var_select <- FSelector::best.first.search(colnames(data[, -1]), eval.fun = function(subset){
      # Computation of the gamma-metric value for the candidate subset of features
      g <- gammaMetric(data[, subset], class = data[, 1], covEstim = estim, plot = FALSE)
      return(g)
    })
    
    return(var_select)
  },
  'GAMMA_FORW' <- function(f, data){
    # Which estimator of the variance-covariance matrix to use
    estim <- ifelse(min(table(data[, 1])) < ncol(data[, -1]), 'shrink', 'empiric')
    
    # Feature selection with the forward search direction
    var_select <- FSelector::forward.search(colnames(data[, -1]), eval.fun = function(subset){
      # Computation of the gamma-metric value for the candidate subset of features
      g <- gammaMetric(data[, subset], class = data[, 1], covEstim = estim, plot = FALSE)
      return(g)
    })
    
    return(var_select)
  },
  'GAMMA_BACK' <- function(f, data){
    # Which estimator of the variance-covariance matrix to use
    estim <- ifelse(min(table(data[, 1])) < ncol(data[, -1]), 'shrink', 'empiric')
    
    # Feature selection with the backward search direction
    var_select <- FSelector::backward.search(colnames(data[, -1]), eval.fun = function(subset){
      # Computation of the gamma-metric value for the candidate subset of features
      g <- gammaMetric(data[, subset], class = data[, 1], covEstim = estim, plot = FALSE)
      return(g)
    })
    
    return(var_select)
  },
  'GAMMA_HC' <- function(f, data){
    # Which estimator of the variance-covariance matrix to use
    estim <- ifelse(min(table(data[, 1])) < ncol(data[, -1]), 'shrink', 'empiric')
    
    # Feature selection with the hill-climbing search direction
    var_select <- FSelector::hill.climbing.search(colnames(data[, -1]), eval.fun = function(subset){
      # Computation of the gamma-metric value for the candidate subset of features
      g <- gammaMetric(data[, subset], class = data[, 1], covEstim = estim, plot = FALSE)
      return(g)
    })
    
    return(var_select)
  }
)


## Simulation iteration -----------------------------------------------------------------------------------------------------------------
# A single function to produce the results of the simulation for a given feature selection method and given train/test data set.
performSimulationIteration <- function(X_train, Y_train, X_test, Y_test, beta, fs){
  ## Inputs: 
  ##    X_train: Matrix of observations X for the train data set (without intercept)
  ##    Y_train: Vectors of labels of the target feature for the train data set
  ##     X_test: Matrix of observations X for the test data set (without intercept)
  ##     Y_test: Vectors of labels of the target feature for the test data set
  ##       beta: Beta vectors without beta0. It was only used to identify quickly informative features
  ##         fs: Name of the feature selection method to apply
  ## Outputs: 
  ##    A dataframe resIteration with multiples column corresponding to the results of one iteration of the simulation on the train and test data set.
  ##      * Approach: Name of the feature selection methods 
  ##      * NbVarSelected: Total number of features selected by the method
  ##      * NbVarInf: Total number of informative features selected 
  ##      * Accuracy_selection: The accuracy of the selection while we defined a True Positive as the feature is informative and is selected by the method
  ##      * Specificity_selection: The specificity of the selection 
  ##      * Sensitivity_selection: The sensitivity of the selection
  ##      * MCC_selection: The MCC of the selection
  ##      * Accuracy_training: The accuracy computed on the adjusted values of the training sample
  ##      * Specificity_training: The specificity computed on the adjusted values of the training sample
  ##      * Sensitivity_training: The sensitivity computed on the adjusted values of the training sample
  ##      * MCC_training: The MCC computed on the adjusted values of the training sample
  ##      * Accuracy_test: The accuracy computed on the predicted values of the test sample
  ##      * Specificity_test: The specificity computed on the predicted values of the test sample
  ##      * Sensitivity_test: The sensitivity computed on the predicted values of the test sample
  ##      * MCC_test: The MCC computed on the predicted values of the test sample
  ##      * VarSelected: The complete list of features selected by the method. To fit everything in a dataframe we pasted all the names of the features
  ##                        in a string separated by ';'
  ##      * coef_estimate: This is the beta coefficients estimated by the logistic regression model. As for VarSelected, we pasted all the values 
  ##                        in a string separated by ';'
  ##      * std_estimate: This is the standard deviation of the estimation of the beta coefficient computed with the logistic regression model. As 
  ##                        for VarSelected, we pasted all the values in a string separated by ';'
  
  # Feature selection step
  var_select <- do.call(fs, list(f = FSelector::as.simple.formula(colnames(X_train), 'y'), data = data.frame(y = factor(Y_train), X_train)))
  
  # Binary vector with 1 being the feature is informative and 0 the feature is non-informative
  non_zero_coefficients <- replace(rep(0, length(beta)), beta != 0, 1)
  
  # Type of features selected by the method
  selected <- (colnames(X_train) %in% var_select)*1
  
  # Computation of the 2x2 confusion matrix for selection
  cm_sel <- caret::confusionMatrix(data = factor(selected, levels = c('0', '1')), reference = factor(non_zero_coefficients), positive = '1')
  
  # Computation of the performance indicators
  spe_sel <- cm_sel$byClass['Specificity']
  sen_sel <- cm_sel$byClass['Sensitivity']
  acc_sel <- cm_sel$overall['Accuracy']
  mcc_sel <- mltools::mcc(preds = factor(selected, levels = c('0', '1')),
                          actuals = factor(non_zero_coefficients))
  
  # Calibration of the Logistic model with the selected features
  ctr <- caret::trainControl(method = 'none', allowParallel = FALSE)
  mod.fit <- caret::train(FSelector::as.simple.formula(var_select, 'y'), 
                          data = data.frame(y = factor(Y_train), X_train), 
                          trControl = ctr, 
                          method = 'glm', 
                          family = 'binomial')
  
  # Get the estimate of the beta coefficients
  resume <- summary(mod.fit)
  coef_estimate <- resume$coef[, 'Estimate']
  std_estimate <- resume$coef[, 'Std. Error']
  
  # Create vector wit NA for non selected features and the estimation of the coefficient beta
  # (or the standard deviation of the estimation) for selected features
  estimates <- data.frame(Beta = rep(NA, length(colnames(X_train))), 
                          Std = rep(NA, length(colnames(X_train))))
  row.names(estimates) <- colnames(X_train)
  estimates[names(coef_estimate)[-1], c('Beta', 'Std')] <- c(coef_estimate[-1], std_estimate[-1])
  
  # Fitted values of the training sample
  fitted <- predict(mod.fit, data.frame(y = Y_train, X_train))
  
  # Classification performances (training sample)
  cm <- caret::confusionMatrix(fitted, factor(Y_train), positive = '1')
  acc_train <- cm$overall['Accuracy']
  spe_train <- cm$byClass['Specificity']
  sen_train <- cm$byClass['Sensitivity']
  mcc_train <- mltools::mcc(fitted, factor(Y_train))
  
  # Prediction of the test sample
  preds <- predict(mod.fit, data.frame(y = Y_test, X_test))
  
  # Classification performances (test sample)
  cm <- caret::confusionMatrix(preds, factor(Y_test), positive = '1')
  acc_test <- cm$overall['Accuracy']
  spe_test <- cm$byClass['Specificity']
  sen_test <- cm$byClass['Sensitivity']
  mcc_test <- mltools::mcc(preds, factor(Y_test))
  
  # Aggregation of the results
  localRes <- data.frame(Approach = fs, NbVarSelected = length(var_select), NbVarInf = sum(non_zero_coefficients & selected), 
                         Accuracy_selection = acc_sel, Specificity_selection = spe_sel, Sensitivity_selection = sen_sel, MCC_selection = mcc_sel, 
                         Accuracy_training = acc_train, Specificity_training = spe_train, Sensitivity_training = sen_train, MCC_training = mcc_train,
                         Accuracy_test = acc_test, Specificity_test = spe_test, Sensitivity_test = sen_test, MCC_test = mcc_test, 
                         VarSelected = paste(var_select, collapse = ';'), coef_estimate = paste0(estimates$Beta, collapse = ';'), std_estimate = paste0(estimates$Std, collapse = ';'))
  return(localRes)
}

## Formating functions ------------------------------------------------------------------------------------------------------------------
dataframeSelectedFeaturesScenario2 <- function(res_table){
  ## Inputs: 
  ##    res_table: The results table of the simulations.
  ## Outputs: 
  ##    VARSELECTED_LONG: The dataframe contains for each method the number of times each feature is selected 
  
  # Number of non-informative features
  p_prime <- res_table$NbVar[1] - n_var_inf
  
  # data.frame with how much time each feature is selected
  VARSELECTED <- data.frame(sapply(Approach, FUN = function(x){
    approach_fs <- which(res$Approach == x)
    list_var_selected <- unlist(strsplit(res$VarSelected[approach_fs], ';'))
    indices <- match(c(paste0('x', 1:n_var_inf), paste0('Noise', 1:n_noise)), names(table(list_var_selected)))
    table(list_var_selected)[indices]
  }))
  
  # The name for SVM-RFE is modified automatically to SVM.RFE in the results, we need to fix this
  colnames(VARSELECTED)[colnames(VARSELECTED) == 'SVM.RFE'] <- 'SVM-RFE'
  
  # Name of the features
  VARSELECTED$Feature <- c(paste0('x', 1:n_var_inf), paste0('N', 1:p_prime))
  
  # Long format
  VARSELECTED_LONG <- reshape2::melt(VARSELECTED, id.vars = 'Feature')
  
  return(VARSELECTED_LONG)
}

dataframeSelectedFeaturesScenario3 <- function(res_table){
  # Function to return how much each feature of Scenario 3 is selected by the feature selection methods
  VARSELECTED <- data.frame(sapply(Approach, FUN = function(x){
    approach_fs <- which(res_table$Approach == x)
    list_var_selected <- unlist(strsplit(res_table$VarSelected[approach_fs], ';'))
    indices <- match(paste0('X', 1:(s_g*n_g)), names(table(list_var_selected)))
    table(list_var_selected)[indices]
  }))
  
  # The name for SVM-RFE is modified automatically to SVM.RFE in the results, we need to change that
  colnames(VARSELECTED)[colnames(VARSELECTED) == 'SVM.RFE'] <- 'SVM-RFE'
  
  # Name of the features
  VARSELECTED$Feature <- paste0('X', 1:(s_g*n_g))
  
  # Long format 
  VARSELECTED_LONG <- reshape2::melt(VARSELECTED, id.vars = 'Feature')
  
  return(VARSELECTED_LONG)
}


