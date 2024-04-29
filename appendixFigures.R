###########################################################################################################################################
####                                                        APPENDICES FIGURES                                                         ####
###########################################################################################################################################

### Code for Figure 8 ---------------------------------------------------------------------------------------------------------------------
source('functionGeneration.R')
set.seed(123)

# Parameters of the generation
n <- 500
n_var_inf <- 2
n_noise <- 1

# Beta values we want to test
beta0 <- c(-2.65, 0, 3, 0, 0, 0)
beta1 <- c(3.6, 3.6, 3.6, 3.6, -1.2, 3.6)
beta2 <- c(-2.2, -2.2, -2.2, -2, -2, 1.2)

# Combinaisons of beta
beta <- data.frame("Combinaison" = LETTERS[1:6], beta0, beta1, beta2)

# Generation of the data with regards to the combinaison of beta values
dat <- do.call('rbind', apply(beta, 1, FUN = function(beta_vec){
  df <- generation(n, n_var_inf, n_noise = 1, c(as.numeric(beta_vec[-1]),0))
  df$Combinaison <- beta_vec[1]
  df$Repartition <- (table(df$y)[2]/n)*100
  df
}))



# Generation of the decision line were rho_i <- 1/2 for each combinaison of beta
decision_line <- do.call('rbind', (apply(beta, 1, FUN = function(beta_vec){
  tmp_beta <- as.numeric(beta_vec[-1])
  x <- seq(-3, 3, 0.01)
  y <- -(tmp_beta[1]/tmp_beta[3]) - (tmp_beta[2]/tmp_beta[3])*x
  df <- data.frame(x = x, y = y, Combinaison = rep(beta_vec[1], length(x)))
  df
})))

# Labeller to get the beta values as the title of each strips
dat$Combinaison <- factor(dat$Combinaison, 
                             levels = LETTERS[1:6], 
                             labels = c(expression(paste(beta[0], " = -2.65 ; ", beta[1], " = 3.6 ; ", beta[2], " = -2.2")),
                                        expression(paste(beta[0], " = 0 ; ", beta[1], " = 3.6 ; ", beta[2], " = -2.2")),
                                        expression(paste(beta[0], " = 3 ; ", beta[1], " = 3.6 ; ", beta[2], " = -2.2")),
                                        expression(paste(beta[0], " = 0 ; ", beta[1], " = 3.6 ; ", beta[2], " = -2")),
                                        expression(paste(beta[0], " = 0 ; ", beta[1], " = -1.2 ; ", beta[2], " = -2")),
                                        expression(paste(beta[0], " = 0 ; ", beta[1], " = 3.6 ; ", beta[2], " = 1.2"))))

# Labeller to get the beta values as the title of each strips
decision_line$Combinaison <- factor(decision_line$Combinaison, 
                                    levels = LETTERS[1:6], 
                                    labels = c(expression(paste(beta[0], " = -2.65 ; ", beta[1], " = 3.6 ; ", beta[2], " = -2.2")),
                                               expression(paste(beta[0], " = 0 ; ", beta[1], " = 3.6 ; ", beta[2], " = -2.2")),
                                               expression(paste(beta[0], " = 3 ; ", beta[1], " = 3.6 ; ", beta[2], " = -2.2")),
                                               expression(paste(beta[0], " = 0 ; ", beta[1], " = 3.6 ; ", beta[2], " = -2")),
                                               expression(paste(beta[0], " = 0 ; ", beta[1], " = -1.2 ; ", beta[2], " = -2")),
                                               expression(paste(beta[0], " = 0 ; ", beta[1], " = 3.6 ; ", beta[2], " = 1.2"))))

# Figure 
plot_figure_08 <- ggplot(data = dat, aes(x = x1, y = x2))+
   geom_point(aes(colour = y), size = 2)+
   geom_text(x = 1.5, y = -2.8, aes(label = paste(Repartition, "% of class 1")), size = 8)+
   geom_line(data = decision_line, aes(x = x, y = y), linewidth = 2)+
   scale_x_continuous(limits = c(-3.1, 3.1))+
   scale_y_continuous(limits = c(-3.1, 3.1))+
   facet_wrap(~Combinaison, ncol = 3, labeller = label_parsed)+
   guides(colour = guide_legend(override.aes = list(size = 4)))+
   theme(strip.text.x = element_text(size = 25), 
         legend.text = element_text(size = 15),
         legend.title = element_text(size = 20))

# Save figure
setEPS()
postscript('Figures/Figure_08.eps', horizontal = FALSE, width = 18, height = 10)
plot_figure_08
dev.off()

### Code for Figure 9 ---------------------------------------------------------------------------------------------------------------------
rm(list = ls())
set.seed(123)

# Let's generate the matrix X with gaussian distribution
n <- 2000
n_var_inf <- 3

# We don't need to generate noisy feature for this example.

# Generation of X with the intercept 
X <- cbind(rep(1, n), 
           matrix(rnorm(n*n_var_inf, mean = 0, sd = 1),
                  ncol = n_var_inf))

# Beta vector for both weak and strong separation (and balanced classes)
beta_weak <- c(0.5, 0.6, -2.5, -1)
beta_strong <- c(0, 3.6, -4, -1)

# Let's compute the rho_i for both cases
rho_i_weak <- exp(X %*% beta_weak)/(1 + exp(X %*% beta_weak))
rho_i_strong <- exp(X %*% beta_strong)/(1 + exp(X %*% beta_strong))
df <- data.frame(Case = rep(c('Weak', 'Strong'), each = n), 
                 rho_i = c(rho_i_weak, rho_i_strong))

# We can now draw the histogram
plot_figure_09 <- ggplot(data = df, aes(x = rho_i))+
  geom_histogram(color = 'grey', fill = 'darkgrey', bins = 50)+
  labs(x = expression(paste(rho[i])))+
  facet_wrap(~Case)+
  theme(strip.text.x = element_text(size = 35),
        axis.text = element_text(size = 30),
        axis.title = element_text(size = 35))

setEPS()
postscript('Figures/Figure_09.eps', horizontal = FALSE, width = 18, height = 10)
plot_figure_09
dev.off()


### Code for Figure 10 --------------------------------------------------------------------------------------------------------------------
rm(list = ls())

# Source wheres the functions are
source('functionGeneration.R')

# Parameters of the covariance matrix  
s_g <- 10                                 # Size of each group
n_g <- 10                                 # Number of groups
constant_cor <- TRUE                      # Constant correlation ?
n_ind <- 1                                # Number of independent group of features

# Parameters of beta
non_zero_coeff <- 1                       # Only the first feature of each group is informative (i.e., is associated with a non-null beta coefficient)
non_zero_group <- c(1, 2, 3, 4, 5, 10)    # First to fith groups have an informative feature and group 10 also.
value <- 1.5                              # Value of the beta for informative features

# Index of informative features
index_of_informative_features <- s_g*(non_zero_group - 1) + non_zero_coeff # Only works because non_zero_coeff is a single value
x_labels <- paste0(1:(s_g*n_g))
x_labels[index_of_informative_features] <- paste0('x', x_labels[index_of_informative_features])
x_labels[-index_of_informative_features] <- paste0('N', x_labels[-index_of_informative_features])

# Construction of the Sigma matrix
rho_1 <- produceRho(s_g = s_g, n_g = n_g, alpha_max = 0.9, n_ind = 1, constant_cor = TRUE)
rho_2 <- produceRho(s_g = s_g, n_g = n_g, alpha_max = 0.6, n_ind = 1, constant_cor = TRUE)
rho_3 <- produceRho(s_g = s_g, n_g = n_g, alpha_max = 0.3, n_ind = 1, constant_cor = TRUE)

colnames(rho_1) = colnames(rho_2) = colnames(rho_3) = x_labels
row.names(rho_1) = row.names(rho_2) = row.names(rho_3) = x_labels

# Figure 

setEPS()
postscript(file = 'Figures/Figure_10.eps', width = 18, height = 10, horizontal = FALSE)
par(mfrow = c(2, 3))
# First plot for alpha_max = 0.9
corrplot::corrplot(rho_1, method = 'square', tl.cex = 0.4, cl.cex = 2, mar = c(0, 0, 2.8, 0), addgrid.col = NA)
title(main = expression(alpha[max] == 0.9), cex.main = 4)

# Second plot for alpha_max = 0.6
corrplot::corrplot(rho_2, method = 'square', tl.cex = 0.4, cl.cex = 2, mar = c(0, 0, 2.8, 0), addgrid.col = NA)
title(main = expression(alpha[max] == 0.6), cex.main = 4)

# Third plot for alpha_max = 0.3
corrplot::corrplot(rho_3, method = 'square', tl.cex = 0.4, cl.cex = 2, mar = c(0, 0, 2.8, 0), addgrid.col = NA)
title(main = expression(alpha[max] == 0.3), cex.main = 4)

# Zoom on the 5th group of features (Sigma_5)
corrplot::corrplot(rho_1[41:50, 41:50], method = 'square', addCoef.col = 'red', number.cex = 2.2, tl.cex = 2, cl.cex = 2)
corrplot::corrplot(rho_2[41:50, 41:50], method = 'square', addCoef.col = 'red', number.cex = 2.2, tl.cex = 2, cl.cex = 2)
corrplot::corrplot(rho_3[41:50, 41:50], method = 'square', addCoef.col = 'red', number.cex = 2.2, tl.cex = 2, cl.cex = 2)

par(mfrow = c(1, 1))
dev.off()

### Code for figure 11 --------------------------------------------------------------------------------------------------------------------
rm(list = ls())

# Non-constant correlation
# Source where the functions are
source('functionGeneration.R')

# Parameters of the covariance matrix  
s_g <- 10                                 # Size of each group
n_g <- 10                                 # Number of groups
constant_cor <- FALSE                     # Constant correlation ?
n_ind <- 1                                # Number of independent group of features

# Parameters of beta
non_zero_coeff <- 1                       # Only the first feature of each group is informative (i.e., is associated with a non-null beta coefficient)
non_zero_group <- c(1, 2, 3, 4, 5, 10)    # First to fith groups have an informative feature and group 10 also.
value <- 1.5                              # Value of the beta for informative features

# Index of informative features
index_of_informative_features <- s_g*(non_zero_group - 1) + non_zero_coeff # Only works because non_zero_coeff is a single value
x_labels <- paste0(1:(s_g*n_g))
x_labels[index_of_informative_features] <- paste0('x', x_labels[index_of_informative_features])
x_labels[-index_of_informative_features] <- paste0('N', x_labels[-index_of_informative_features])

# Construction of the Sigma matrix for 3 different levels of maximum correlation alpha_max
rho_1 <- produceRho(s_g = s_g, n_g = n_g, alpha_max = 0.9, c = 0.35, n_ind = n_ind, constant_cor = constant_cor)
rho_2 <- produceRho(s_g = s_g, n_g = n_g, alpha_max = 0.6, c = 0.25, n_ind = n_ind, constant_cor = constant_cor)
rho_3 <- produceRho(s_g = s_g, n_g = n_g, alpha_max = 0.3, c = 0.1, n_ind = n_ind, constant_cor = constant_cor)

colnames(rho_1) = colnames(rho_2) = colnames3 = x_labels
row.names(rho_1) = row.names(rho_2) = row.names(rho_3) = x_labels

# Figure
setEPS()
postscript('Figures/Figure_11.eps', width = 18, height = 10, horizontal = FALSE)
par(mfrow = c(2, 3))

# First figure for high correlation alpha_max = 0.9 and c = 0.35
corrplot::corrplot(rho_1, method = 'square', tl.cex = 0.4, cl.cex = 2, mar = c(0,0,2.8,0), addgrid.col = NA)
title(main = expression(paste(alpha[max] == 0.9, ' and ', c == 0.35)), cex.main = 3)

# Second figure for medium correlation alpha_max = 0.6 and c = 0.25
corrplot::corrplot(rho_2, method = 'square', tl.cex = 0.4, cl.cex = 2, mar = c(0,0,2.8,0), addgrid.col = NA)
title(main = expression(paste(alpha[max] == 0.6, ' and ', c == 0.25)), cex.main = 3)

# Third figure for low correlation alpha_max = 0.3 and c = 0.1
corrplot::corrplot(rho_3, method = 'square', tl.cex = 0.4, cl.cex = 2, mar = c(0,0,2.8,0), addgrid.col = NA)
title(main = expression(paste(alpha[max] == 0.3, ' and ', c == 0.1)), cex.main = 3)

# Zoom on the 5th group of features
corrplot::corrplot(rho_1[41:50, 41:50], method = 'square', addCoef.col = 'red', number.cex = 1.8, tl.cex = 2, cl.cex = 2)
corrplot::corrplot(rho_2[41:50, 41:50], method = 'square', addCoef.col = 'red', number.cex = 1.8, tl.cex = 2, cl.cex = 2)
corrplot::corrplot(rho_3[41:50, 41:50], method = 'square', addCoef.col = 'red', number.cex = 1.8, tl.cex = 2, cl.cex = 2)

par(mfrow = c(1, 1))
dev.off()

### Code for Figure 12 --------------------------------------------------------------------------------------------------------------------------
rm(list = ls())

# Source of the functions for generation of the data
source('functionGeneration.R')

# Parameters of the generation of beta vector
s_g <- 10                                    # Size of each group of features
n_g <- 10                                    # Number of group of features
non_zero_coeff <- 1                          # Index of informative features in each group
non_zero_group <- c(1, 2, 3, 4, 5, 10)       # Index of the group with informative features
value <- 1.5                                 # Value of the beta for informative features

# Generation of the vector
tmp_res <- produceBeta(s_g = s_g, n_g = n_g, non_zero_coeff = non_zero_coeff, non_zero_group = non_zero_group, value = value)
beta0 <- tmp_res[1]
beta <- tmp_res[-1]
betaDataFrame <- data.frame(Beta = beta, Group = rep(1:n_g, each = s_g), Index = rep(1:s_g, times = n_g))

plot_figure_12 <- ggplot(data = betaDataFrame, aes(factor(Index), rev(Group)))+
  geom_tile(aes(fill = factor(Beta)), color = 'white', lwd = 1.5, linetype = 1)+
  geom_text(aes(label = paste('beta', '[', 1:(s_g*n_g), ']')), parse = TRUE, size = 15)+
  scale_y_discrete(labels = rev(1:10), limits = factor(rev(1:10)))+
  labs(fill = expression(paste(beta, ' value')), x = 'Index', y = 'Group')+
  theme_minimal()+
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 30),
        legend.key.size = unit(1, units = 'cm'), 
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30))

setEPS()
postscript('Figures/Figure_12.eps', width = 18, height = 10, horizontal = FALSE)
plot_figure_12
dev.off()




###########################################################################################################################################
####                                                   SUPPORTING INFORMATION                                                          ####
###########################################################################################################################################
# Uploading the results of scenario 2
rm(list = ls())

# Function
source('functions.R')

## Parameters of the generation 
N <- 2000                                          # Number of observations
n_var_inf <- 3                                     # Number of informative features
n_noise <- c(47, 97, 147, 197)                     # Number of non-informative features
separation <- c('faible', 'forte')                 # Separation between classes
repartition <- c('equilibre', 'desequilibre')      # Repartition between classes
R <- 50                                             # Number of repetition

# Feature selection methods
Approach <- c('BASELINE', 'BEST', 'CFS', 'CHI2', 'CONS', 'IG', 'IGR', 'ONER', 
              'RELIEF', 'RFI', 'SU', 'SVM-RFE', 'GAMMA_BACK', 'GAMMA_BF', 'GAMMA_FORW', 'GAMMA_HC')

# Upload the results
# All the files in scenario 2
files <- dir('Scenario2/')

# Get files related to the same number of features 
files_50 <- paste('Scenario2/', files[grep('res_50_repetitions_50', files)], sep = '')
files_100 <- paste('Scenario2/', files[grep('res_50_repetitions_100', files)], sep = '')
files_150 <- paste('Scenario2/', files[grep('res_50_repetitions_150', files)], sep = '')

# Aggregating the files in one dataframe
res_50 <- data.frame(do.call('rbind', lapply(files_50, read.table)))
res_100 <- data.frame(do.call('rbind', lapply(files_100, read.table)))
res_150 <- data.frame(do.call('rbind', lapply(files_150, read.table)))

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
  
}

#### 50 features ----------------------------------------------------------------------------------------------------------------------------------------------------------------
## Code for Table S1 --------------------------------------------------------------------------------------------------------------------------------------------------------

# Upper part
print(xtable::xtable(cbind(RES_50[which(RES_50$Separation == 'forte' & RES_50$Repartition == 'equilibre'), 
             c('Approach', 'NbVarSelected', 'NbVarInf', 'Specificity_selection', 'Sensitivity_selection', 'MCC_test')],
      RES_50[which(RES_50$Separation == 'faible' & RES_50$Repartition == 'equilibre'),
             c('NbVarSelected', 'NbVarInf', 'Specificity_selection', 'Sensitivity_selection', 'MCC_test')]), style = 'latex', digits = c(0, 0, 2, 2, 2, 2, 3, 2, 2, 2, 2, 3)))

# Lower part
print(xtable::xtable(cbind(RES_50[which(RES_50$Separation == 'forte' & RES_50$Repartition == 'desequilibre'),
             c('Approach', 'NbVarSelected', 'NbVarInf', 'Specificity_selection', 'Sensitivity_selection', 'MCC_test')],
      RES_50[which(RES_50$Separation == 'faible' & RES_50$Repartition == 'desequilibre'),
             c('NbVarSelected', 'NbVarInf', 'Specificity_selection', 'Sensitivity_selection', 'MCC_test')]), style = 'latex', digits = c(0, 0, 2, 2, 2, 2, 3, 2, 2, 2, 2, 3)))

## Code for Figure S1 ------------------------------------------------------------------------------------------------------------------------------------------------------
# Dataframe of selected features for each situation with 50 features
VARSELECTED_LONG_faible_desequilibre <- dataframeSelectedFeaturesScenario2(res_50[which(res_50$Separation == 'faible' & res_50$Repartition == 'desequilibre'),])
VARSELECTED_LONG_faible_equilibre <- dataframeSelectedFeaturesScenario2(res_50[which(res_50$Separation == 'faible' & res_50$Repartition == 'equilibre'),])
VARSELECTED_LONG_forte_desequilibre <- dataframeSelectedFeaturesScenario2(res_50[which(res_50$Separation == 'forte' & res_50$Repartition == 'desequilibre'),])
VARSELECTED_LONG_forte_equilibre <- dataframeSelectedFeaturesScenario2(res_50[which(res_50$Separation == 'forte' & res_50$Repartition == 'equilibre'),])

# Modification of the labels
VARSELECTED_LONG_faible_desequilibre[, c('Balance', 'Separation')] <- rep(c('Unbalanced', 'Weak'), each = nrow(VARSELECTED_LONG_faible_desequilibre))
VARSELECTED_LONG_faible_equilibre[, c('Balance', 'Separation')] <- rep(c('Balanced', 'Weak'), each = nrow(VARSELECTED_LONG_faible_equilibre))
VARSELECTED_LONG_forte_desequilibre[, c('Balance', 'Separation')] <- rep(c('Unbalanced', 'Strong'), each = nrow(VARSELECTED_LONG_forte_desequilibre))
VARSELECTED_LONG_forte_equilibre[, c('Balance', 'Separation')] <- rep(c('Balanced', 'Strong'), each = nrow(VARSELECTED_LONG_forte_equilibre))

# Binding to a single dataframe
VARSELECTED_LONG <- rbind(VARSELECTED_LONG_faible_desequilibre,
                          VARSELECTED_LONG_faible_equilibre,
                          VARSELECTED_LONG_forte_desequilibre,
                          VARSELECTED_LONG_forte_equilibre)

# Plot 
plot_figure_S1 <- ggplot(data = VARSELECTED_LONG, aes(x = Feature, y = variable))+
  geom_tile(aes(fill = value))+
  scale_x_discrete(breaks = unique(VARSELECTED_LONG$Feature), limits = unique(VARSELECTED_LONG$Feature))+
  scale_y_discrete(limits = rev(Approach))+
  scale_fill_gradientn(colours = rev(viridis::inferno(10)), breaks = seq(0, R, R/5))+
  labs(fill = '', x = '', y = '')+
  facet_grid(Balance~Separation)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 14, hjust = 1),
        axis.text.y = element_text(size = 20),
        legend.key.width = unit(0.75, 'cm'), 
        legend.key.height = unit(2, 'cm'),
        legend.margin = margin(0.1, 0.1, 0.1, 0.1),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 30, margin = margin(0.025, 0, 0.025, 0, 'cm')),
        strip.text.y = element_text(size = 30, margin = margin(0, 0.025, 0, 0.025, 'cm')))

# Save Figure
setEPS()
postscript('Figures/Figure_S1.eps', width = 18, height = 10, horizontal = FALSE)
plot_figure_S1
dev.off()

## Code for Figure S2 ------------------------------------------------------------------------------------------------------------------------------------------------------
# Changing the labels 
res_50$Separation <- factor(res_50$Separation, levels = c('forte', 'faible'), labels = c('Strong', 'Weak'))
res_50$Repartition <- factor(res_50$Repartition, levels = c('equilibre', 'desequilibre'), labels = c('Balanced', 'Unbalanced'))

# Plot 
plot_figure_S2 <- ggplot(data = res_50, aes(x = Approach, y = MCC_test))+
  geom_boxplot(size = 0.3, outlier.size = 0.2, outlier.shape = 19)+
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

# Save
setEPS()
postscript('Figures/Figure_S2.eps', width = 18, height = 10, horizontal = FALSE)
plot_figure_S2
dev.off()


#### 100 features ----------------------------------------------------------------------------------------------------------------------------------------------------------------

## Code for Table S2 -------------------------------------------------------------------------------------------------------------------------------------------------------------

# Upper part
print(xtable::xtable(cbind(RES_100[which(RES_100$Separation == 'forte' & RES_100$Repartition == 'equilibre'),
              c('Approach', 'NbVarSelected', 'NbVarInf', 'Specificity_selection', 'Sensitivity_selection', 'MCC_test')],
      RES_100[which(RES_100$Separation == 'faible' & RES_100$Repartition == 'equilibre'),
              c('NbVarSelected', 'NbVarInf', 'Specificity_selection', 'Sensitivity_selection', 'MCC_test')]), style = 'latex', digits = c(0, 0, 2, 2, 2, 2, 3, 2, 2, 2, 2, 3)))

# Lower part
print(xtable::xtable(cbind(RES_100[which(RES_100$Separation == 'forte' & RES_100$Repartition == 'desequilibre'), 
              c('Approach', 'NbVarSelected', 'NbVarInf', 'Specificity_selection', 'Sensitivity_selection', 'MCC_test')],
      RES_100[which(RES_100$Separation == 'faible' & RES_100$Repartition == 'desequilibre'), 
              c('NbVarSelected', 'NbVarInf', 'Specificity_selection', 'Sensitivity_selection', 'MCC_test')]), style = 'latex', digits = c(0, 0, 2, 2, 2, 2, 3, 2, 2, 2, 2, 3)))

## Code for figure S3 ------------------------------------------------------------------------------------------------------------------------------------------------------------

VARSELECTED_LONG_faible_desequilibre <- dataframeSelectedFeaturesScenario2(res_100[which(res_100$Separation == 'faible' & res_100$Repartition == 'desequilibre'),])
VARSELECTED_LONG_faible_equilibre <- dataframeSelectedFeaturesScenario2(res_100[which(res_100$Separation == 'faible' & res_100$Repartition == 'equilibre'),])
VARSELECTED_LONG_forte_desequilibre <- dataframeSelectedFeaturesScenario2(res_100[which(res_100$Separation == 'forte' & res_100$Repartition == 'desequilibre'),])
VARSELECTED_LONG_forte_equilibre <- dataframeSelectedFeaturesScenario2(res_100[which(res_100$Separation == 'forte' & res_100$Repartition == 'equilibre'),])

VARSELECTED_LONG_faible_desequilibre[, c('Balance', 'Separation')] <- rep(c('Unbalanced', 'Weak'), each = nrow(VARSELECTED_LONG_faible_desequilibre))
VARSELECTED_LONG_faible_equilibre[, c('Balance', 'Separation')] <- rep(c('Balanced', 'Weak'), each = nrow(VARSELECTED_LONG_faible_equilibre))
VARSELECTED_LONG_forte_desequilibre[, c('Balance', 'Separation')] <- rep(c('Unbalanced', 'Strong'), each = nrow(VARSELECTED_LONG_forte_desequilibre))
VARSELECTED_LONG_forte_equilibre[, c('Balance', 'Separation')] <- rep(c('Balanced', 'Strong'), each = nrow(VARSELECTED_LONG_forte_equilibre))

VARSELECTED_LONG <- rbind(VARSELECTED_LONG_faible_desequilibre, VARSELECTED_LONG_faible_equilibre, 
                          VARSELECTED_LONG_forte_desequilibre, VARSELECTED_LONG_forte_equilibre)

plot_figure_S3 <- ggplot(data = VARSELECTED_LONG, aes(x = Feature, y = variable))+
  geom_tile(aes(fill = value))+
  scale_x_discrete(breaks = unique(VARSELECTED_LONG$Feature), limits = unique(VARSELECTED_LONG$Feature))+
  scale_fill_gradientn(colours = rev(viridis::inferno(10)), breaks = seq(0, R, R/5))+
  scale_y_discrete(limits = rev(Approach))+
  labs(fill = '', x = '', y = '')+
  facet_grid(Balance~Separation)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5, hjust = 1),
        axis.text.y = element_text(size = 20),
        legend.key.width = unit(0.75, 'cm'),
        legend.key.height = unit(2, 'cm'),
        legend.margin = margin(0.1, 0.1, 0.1, 0.1), 
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 25, margin = margin(0.025, 0, 0.025, 0, "cm")),
        strip.text.y = element_text(size = 25, margin = margin(0, 0.025, 0, 0.025, 'cm')))

setEPS()
postscript('Figures/Figure_S3.eps', width = 18, height = 10, horizontal = FALSE)
plot_figure_S3
dev.off()

## Code for figure S4 ------------------------------------------------------------------------------------------------------------------------------------------------------------
res_100$Separation <- factor(res_100$Separation, levels = c('forte', 'faible'), labels = c('Strong', 'Weak'))
res_100$Repartition <- factor(res_100$Repartition, levels = c('equilibre', 'desequilibre'), labels = c('Balanced', 'Unbalanced'))

plot_figure_S4 <- ggplot(data = res_100, aes(x = Approach, y = MCC_test))+
  geom_boxplot(size = 0.3, outlier.size = 0.2, outlier.shape = 19)+
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
postscript('Figures/Figure_S4.eps', width = 18, height = 10, horizontal = FALSE)
plot_figure_S4
dev.off()

#### 150 features -------------------------------------------------------------------------------------------------------------------------------------------------------------
## Code for Table S3 ----------------------------------------------------------------------------------------------------------------------------------------------------------

# Upper part
print(xtable::xtable(cbind(RES_150[which(RES_150$Separation == 'forte' & RES_150$Repartition == 'equilibre'),
              c('Approach', 'NbVarSelected', 'NbVarInf', 'Specificity_selection', 'Sensitivity_selection', 'MCC_test')],
      RES_150[which(RES_150$Separation == 'faible' & RES_150$Repartition == 'equilibre'),
              c('NbVarSelected', 'NbVarInf', 'Specificity_selection', 'Sensitivity_selection', 'MCC_test')]), style = 'latex', digits = c(0, 0, 2, 2, 2, 2, 3, 2, 2, 2, 2, 3)))

# Lower part
print(xtable::xtable(cbind(RES_150[which(RES_150$Separation == 'forte' & RES_150$Repartition == 'desequilibre'),
              c('Approach', 'NbVarSelected', 'NbVarInf', 'Specificity_selection', 'Sensitivity_selection', 'MCC_test')],
      RES_150[which(RES_150$Separation == 'faible' & RES_150$Repartition == 'desequilibre'),
              c('NbVarSelected', 'NbVarInf', 'Specificity_selection', 'Sensitivity_selection', 'MCC_test')]), style = 'latex', digits = c(0, 0, 2, 2, 2, 2, 3, 2, 2, 2, 2, 3)))

## Code for figure S5 ---------------------------------------------------------------------------------------------------------------------------------------------------------
VARSELECTED_LONG_faible_desequilibre <- dataframeSelectedFeaturesScenario2(res_150[which(res_150$Separation == 'faible' & res_150$Repartition == 'desequilibre'),])
VARSELECTED_LONG_faible_equilibre <- dataframeSelectedFeaturesScenario2(res_150[which(res_150$Separation == 'faible' & res_150$Repartition == 'equilibre'),])
VARSELECTED_LONG_forte_desequilibre <- dataframeSelectedFeaturesScenario2(res_150[which(res_150$Separation == 'forte' & res_150$Repartition == 'desequilibre'),])
VARSELECTED_LONG_forte_equilibre <- dataframeSelectedFeaturesScenario2(res_150[which(res_150$Separation == 'forte' & res_150$Repartition == 'equilibre'),])

VARSELECTED_LONG_faible_desequilibre[, c('Balance', 'Separation')] <- rep(c('Unbalanced', 'Weak'), each = nrow(VARSELECTED_LONG_faible_desequilibre))
VARSELECTED_LONG_faible_equilibre[, c('Balance', 'Separation')] <- rep(c('Balanced', 'Weak'), each = nrow(VARSELECTED_LONG_faible_equilibre))
VARSELECTED_LONG_forte_desequilibre[, c('Balance', 'Separation')] <- rep(c('Unbalanced', 'Strong'), each = nrow(VARSELECTED_LONG_forte_desequilibre))
VARSELECTED_LONG_forte_equilibre[, c('Balance', 'Separation')] <- rep(c('Balanced', 'Strong'), each = nrow(VARSELECTED_LONG_forte_equilibre))

VARSELECTED_LONG <- rbind(VARSELECTED_LONG_faible_desequilibre, VARSELECTED_LONG_faible_equilibre, 
                          VARSELECTED_LONG_forte_desequilibre, VARSELECTED_LONG_forte_equilibre)

plot_figure_S5 <- ggplot(data = VARSELECTED_LONG, aes(x = Feature, y = variable))+
  geom_tile(aes(fill = value))+
  scale_x_discrete(breaks = unique(VARSELECTED_LONG$Feature), limits = unique(VARSELECTED_LONG$Feature))+
  scale_fill_gradientn(colours = rev(viridis::inferno(10)), breaks = seq(0, R, R/5))+
  scale_y_discrete(limits = rev(Approach))+
  labs(fill = '', x = '', y = '')+
  facet_grid(Balance~Separation)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 3.5, hjust = 1),
        axis.text.y = element_text(size = 20),
        legend.key.width = unit(0.75, 'cm'),
        legend.key.height = unit(2, 'cm'),
        legend.margin = margin(0.1, 0.1, 0.1, 0.1), 
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 25, margin = margin(0.025, 0, 0.025, 0, "cm")),
        strip.text.y = element_text(size = 25, margin = margin(0, 0.025, 0, 0.025, 'cm')))

setEPS()
postscript('Figures/Figure_S5.eps', width = 18, height = 10, horizontal = FALSE)
plot_figure_S5
dev.off()

## Code for Figure S6 -------------------------------------------------------------------------------------------------------------------------------------------------------
res_150$Separation <- factor(res_150$Separation, levels = c('forte', 'faible'), labels = c('Strong', 'Weak'))
res_150$Repartition <- factor(res_150$Repartition, levels = c('equilibre', 'desequilibre'), labels = c('Balanced', 'Unbalanced'))

plot_figure_S6 <- ggplot(data = res_150, aes(x = Approach, y = MCC_test))+
  geom_boxplot(size = 0.3, outlier.size = 0.2, outlier.shape = 19)+
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
postscript('Figures/Figure_S6.eps', width = 18, height = 10, horizontal = FALSE)
plot_figure_S6
dev.off()






