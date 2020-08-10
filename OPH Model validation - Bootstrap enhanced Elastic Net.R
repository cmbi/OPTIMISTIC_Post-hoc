# Script run with R version 3.5.1 (2018-07-12)
# Last changes applied on: 30/07/2020

# Script used to validate the Bootstrap enhanced Elsatic-Net approach for variable selection in the OPTIMISTIC dataset

set.seed(1)

####################
## Load libraries ##
####################

library("readxl")       #V 1.3.1
library("dplyr")        #V 0.8.5
library("MASS")         #V 7.3.51.5
library("glmnet")       #V 2.0.16
library("boot")         #V 1.3.20
library("ggplot2")      #V 3.1.1
library("gridExtra")    #V 2.3
library("reshape2")     #V 1.4.4


###############
## Functions ##
###############

## Function to run cv.glmnet on bootstrap distributions and return coefficients
fitglmnet <- function(data, indices, alpha, nfolds){
  # Allow boot function to select sample
  d <- data[indices,]
  
  # Matrix conversion & splitting dependent variable
  x <- as.matrix(d[,-which(names(d)=="DM1Activc")])
  y <- as.matrix(d[,"DM1Activc", drop=F])
  
  # Fit model and return coefficients
  fit <- cv.glmnet(x, y, type.measure ="mse", alpha = alpha, family="gaussian", nfolds = nfolds)
  coef <- coef(fit, s="lambda.1se")
  return(as.vector(coef))
}

## Function to build a linear model based on bootstrap predictor frequency distributions
lm_bs_enet <- function(boot_freq, R, freq, data){
  
  # Calculate variable inclusion frequencies
  boot_sums <- colSums(boot_freq != 0) / R
  names(boot_sums) <- gsub(x = names(boot_sums), pattern = ":", replacement = "*")
  
  # Build linear model based on selected variables
  boot_predictors <- names(boot_sums[boot_sums >= freq])
  boot_predictors <- boot_predictors[boot_predictors != "Intercept"]
  boot_predictors <- paste(boot_predictors, collapse = " + ")
  boot_formula <- paste("DM1Activc", boot_predictors, sep=" ~ ")
  lm_boot <- lm(as.formula(boot_formula), data=data)
  return(lm_boot)
}

## Function to calculate r^2 and adjusted r^2 on testingset / prediction data
# p = number of predictors excluding constant; N = total sample size
rsq <- function(observed, predicted, p, N){
  # r^2
  SSTot <- sum((observed - mean(observed))^2)
  SSRes <- sum((observed - predicted)^2)
  Rsq <- round(1 - (SSRes / SSTot),2)
  
  # Adjusted r^2
  Rsq_adj <- round(1 - (((1-Rsq) * (N - 1)) / (N - p - 1)),2)
  return(c(Rsq, Rsq_adj))
}


############################################
## Set working directory and open databse ##
############################################
setwd("Z:/documents")
df <- read_excel("Outcome measures overview (V5)_corrected.xlsx", sheet = 1)


###############################################
## Create dataframe with relevant predictors ##
###############################################

## Subset relevant predictors for somatic instability calculation & convert variables to numeric
pdf <- data.frame(df$LogV2Diff, df$LogePAL, df$VariantRepeats, Age=df$AgeV2)
colnames(pdf) <- gsub(x = colnames(pdf), pattern = "df.", replacement = "")
pdf <- sapply(pdf[,1:length(pdf)], function(x){as.numeric(as.character(x))})
pdf <- as.data.frame(pdf)

## Add patient ID's & select complete cases | n=247, 8 cases removed
pdf$PatientID <- df$PatientID...164
pdf <- pdf[complete.cases(pdf),]  

## Calculate somatic instability (SI) & add to prediction dataframe
SI <- lm(LogV2Diff ~ LogePAL + Age + LogePAL*Age + VariantRepeats, data=pdf)
pdf$SI <- rstandard(SI)

## Adding baseline DM1-Activ-c, Sex and treatment to pdf and select complete cases | n=246, 1 case removed
subset <- data.frame(PatientID=df$PatientID...164, DM1Activc=df$DM1ActivCV2, Sex=df$SexCode, treatment=df$TreatmentCode)
pdf <- merge(pdf, subset, by="PatientID")
pdf <- pdf[complete.cases(pdf),]  

## Remove Patient ID's, LogV2DIFF and treatment
# although removed as a predictor, treatment status is retained to allow for later subsetting
pdf$PatientID <- NULL
pdf$LogV2Diff <- NULL
treatment <- pdf$treatment
pdf$treatment <- NULL

## Calculation of interaction variables
pdf$'Age:LogePAL' <- pdf$LogePAL*pdf$Age
pdf$'LogePAL:SI' <- pdf$LogePAL*pdf$SI
pdf$'Age:SI' <- pdf$SI*pdf$Age
pdf$'Age:LogePAL:SI' <- pdf$LogePAL*pdf$SI*pdf$Age


################################
## Implementation of analyses ##
################################

# Dataframe to store selected predictors
selected_pred <- data.frame(matrix(ncol=10, nrow=0))
colnames(selected_pred) <- c("(Intercept)", colnames(pdf)[-5])

# Dataframe to store mean statistical estimates
mean_est <- data.frame(matrix(ncol=4, nrow=0))
colnames(mean_est) <- c("IS.a.Rsq", "IS.RMSE", "OS.a.Rsq", "OS.RMSE")


for (x in c("Full", "Intervention")){
  
  ## Subset Intervention cohort in 2nd loop
  if (x=="Intervention"){
    pdf <- pdf[treatment == 1,]
    plottags <- c("B", "D")
  } else {
    plottags <- c("A", "C")
  }
  
  ####################################
  ## Backwards stepwise regression  ##
  ####################################
  
  ## Build the model
  lm <- lm(DM1Activc ~ Age + LogePAL + SI + VariantRepeats + Sex + 
             LogePAL*Age + LogePAL*SI + SI*Age + LogePAL*SI*Age, data=pdf)
  lm_step <- stepAIC(lm, direction = "backward")
  summary(lm_step) #exact reproduction of published data
  
  ## Store selected variables
  # Selected coefficients
  step_pred <- lm_step$coefficients
  step_pred[!is.na(step_pred)] <- 1
  selected_pred <- bind_rows(selected_pred, step_pred)
  
  # Statistical estimates
  lm_step_radj <- round(summary(lm_step)$adj.r.squared,2)
  lm_step_rmse <- round(sqrt(mean(lm_step$residuals^2)),2)
  temp <- data.frame('IS.a.Rsq'=lm_step_radj, 'IS.RMSE'=lm_step_rmse)
  mean_est <- bind_rows(mean_est, temp)
  
  
  ##################################################################
  ## Bootstrap enhanced Elastic-Net regression variable selection ##
  ##################################################################
  
  ## Fit Elastic-Net on 5000 bootstrap distributions | 5 fold CV, alpha=0.5 
  boot_res <- boot(data=pdf, fitglmnet, R=5000, alpha=0.5, nfolds=5)
  boot_freq <- as.data.frame(boot_res$t)
  colnames(boot_freq) <- c("Intercept", colnames(pdf)[-5])  
  
  ## Build linear model based on predictor inclusion frequencies
  lm_boot <- lm_bs_enet(boot_freq, 5000, 0.6, pdf)
  summary(lm_boot)
  
  ## Store selected variables
  # Selected coefficients
  boot_pred <- lm_boot$coefficients
  boot_pred[!is.na(boot_pred)] <- 1
  names(boot_pred) <- gsub(x = names(boot_pred), pattern = "LogePAL:Age", replacement = "Age:LogePAL")
  selected_pred <- bind_rows(selected_pred, boot_pred)
  
  # Statistical estimates
  lm_boot_radj <- round(summary(lm_boot)$adj.r.squared,2)
  lm_boot_rmse <- round(sqrt(mean(lm_boot$residuals^2)),2)
  temp <- data.frame('IS.a.Rsq'=lm_boot_radj, 'IS.RMSE'=lm_boot_rmse)
  mean_est <- bind_rows(mean_est, temp)
  
  
  #################################################
  ## Model validation on training/testing splits ##
  #################################################
  step_coefficients <- data.frame() 
  boot_coefficients <- data.frame()
  radj_estimates <- data.frame()
  rsme_estimates <- data.frame()
  
  
  for (i in 1:10){
    #####################################################################################
    ## Randomly splitting dataset into 75% training (n=184|94) / 25% testing (n=62|32) ##
    #####################################################################################
    rrows <- sample(1:nrow(pdf), 0.75*nrow(pdf))
    train <- pdf[rrows,]
    test <- pdf[-rrows,]
    
    ###################################
    ## Backwards stepwise regression ##
    ###################################
    
    ## Build model using stepAIC function
    lm <- lm(DM1Activc ~ Age + LogePAL + SI + VariantRepeats + Sex + 
               LogePAL*Age + LogePAL*SI + SI*Age + LogePAL*SI*Age, data=train)
    lm_step <- stepAIC(lm, direction = "backward")
    
    ## Store selected coefficients
    step_coefficients <- bind_rows(step_coefficients, lm_step$coefficients)
    
    ## Obtain statistical estimates on training set
    train_step_radj <- round(summary(lm_step)$adj.r.squared,2)
    train_step_rmse <- round(sqrt(mean(lm_step$residuals^2)),2)
    
    ## Predict testing set and obtain statistical estimates
    p = (length(lm_step$coefficients) - 1) # number of predictors for radj calculation; excluding intercept
    test$lm_step_predict <- predict(lm_step, newdata = test)
    test_step_radj <- rsq(test$DM1Activc, test$lm_step_predict, p, nrow(test))[2]
    test_step_rmse <- round(sqrt(mean((test$DM1Activc - test$lm_step_predict)^2)),2)
    
    
    ###############################################
    ## Bootstrap enhanced Elastic-Net regression ##
    ###############################################
    
    ## Fit Elastic-Net on 5000 bootstrap distributions | 5 fold CV, alpha=0.5  
    boot_res <- boot(data=pdf, fitglmnet, R=5000, alpha=0.5, nfolds=5)
    boot_freq <- as.data.frame(boot_res$t)
    colnames(boot_freq) <- c("Intercept", colnames(pdf)[-5])
    
    ## Build linear model based on predictor inclusion frequencies
    lm_boot <- lm_bs_enet(boot_freq, 5000, 0.6, train)
    
    ## Store selected coefficients
    boot_pred <- lm_boot$coefficients
    names(boot_pred) <- gsub(x = names(boot_pred), pattern = "LogePAL:Age", replacement = "Age:LogePAL")
    boot_coefficients <- bind_rows(boot_coefficients, boot_pred)
    
    ## Obtain statistical estimates on training set
    train_boot_radj <- round(summary(lm_boot)$adj.r.squared,2)
    train_boot_rmse <- round(sqrt(mean(lm_boot$residuals^2)),2)
    
    ## Predict testing set and obtain statistical estimates
    p = (length(lm_boot$coefficients) - 1)
    test$lm_boot_predict <- predict(lm_boot, newdata = test)
    test_boot_radj <- rsq(test$DM1Activc, test$lm_boot_predict, p, nrow(test))[2]
    test_boot_rmse <- round(sqrt(mean((test$DM1Activc - test$lm_boot_predict)^2)),2)
    
    
    #############################################
    ## Combine and store statistical estimates ##
    #############################################
    
    radj_temp <- data.frame(train_step_radj, train_boot_radj, test_step_radj, test_boot_radj)
    radj_estimates <- rbind(radj_estimates, radj_temp)
    
    rmse_temp <- data.frame(train_step_rmse, train_boot_rmse, test_step_rmse, test_boot_rmse)
    rsme_estimates <- rbind(rsme_estimates, rmse_temp)
    
  }
  
  ############################################################################
  ## Calculate variable inclusion frequncies and add to predictor dataframe ##
  ############################################################################
  
  step_sums <- colSums(!is.na(step_coefficients))
  step_sums <- round(step_sums/10, 2)
  selected_pred <- bind_rows(selected_pred, step_sums)
  
  boot_sums <- colSums(!is.na(boot_coefficients))
  boot_sums <- round(boot_sums/10, 2)
  selected_pred <- bind_rows(selected_pred, boot_sums)
  selected_pred[is.na(selected_pred)] <- 0
  
  
  ###################################################
  ## Visualization of a-rsq and rmse distributions ##
  ###################################################

  ## Update variable names for visualization
  names(radj_estimates) <- c("BSR-IS", "BeEN-IS", "BSR-OS", "BeEN-OS")
  names(rsme_estimates) <- c("BSR-IS", "BeEN-IS", "BSR-OS", "BeEN-OS")
  
  assign(paste(x, "radj", sep="_"),   
         ggplot(data= melt(radj_estimates), aes(x=variable, y=value))+
           geom_boxplot() +
           scale_y_continuous(expand = c(0,0), limits = c(-1.6, 0.3), breaks = scales::pretty_breaks(n=10)) +
           ggtitle(paste(x, "cohort")) +
           labs(x="", y="Adjusted R-squared", tag=plottags[1]) +
           theme(
             panel.border = element_rect(colour="black", fill = NA, size = 0.4),
             axis.text = element_text(size = 15),
             axis.title = element_text(colour="black", size = 17, face="bold"),
             plot.title = element_text(colour="black", size = 17, face="bold"),
             plot.tag.position = "topleft",
             plot.tag = element_text(size=25, face="bold"),
             panel.background = element_rect(fill="white"),
             aspect.ratio = 1,
             plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm")
           ))
  
  assign(paste(x, "rsme", sep="_"), 
         ggplot(data= melt(rsme_estimates), aes(x=variable, y=value))+
           geom_boxplot() +
           labs(x="", y="Root Mean Square Error", tag=plottags[2]) +
           scale_y_continuous(expand = c(0,0), limits = c(10, 21), breaks = scales::pretty_breaks(n=10)) +
           ggtitle(paste(x, "cohort")) +
           theme(
             panel.border = element_rect(colour="black", fill = NA, size = 0.4),
             axis.text = element_text(size = 15),
             axis.title = element_text(colour="black", size = 17, face="bold"),
             plot.title = element_text(colour="black", size = 17, face="bold"),
             plot.tag.position = "topleft",
             plot.tag = element_text(size=25, face="bold"),
             panel.background = element_rect(fill="white"),
             aspect.ratio = 1,
             plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm")
           ))
  
  ############################################################
  ## Calculate mean a-rsq and rsme values and store results ## 
  ############################################################
  
  radj_estimates <- rbind(radj_estimates, colMeans(radj_estimates))
  row.names(radj_estimates) <- c(1:10, "Mean")
  
  rsme_estimates <- rbind(rsme_estimates, colMeans(rsme_estimates))
  row.names(rsme_estimates) <- c(1:10, "Mean")
  
  assign(paste(x,"radj_estimates",sep = "_"), radj_estimates)
  assign(paste(x,"rsme_estimates",sep = "_"), rsme_estimates)
  
  ##############################################
  ## Add mean estimates to mean_est dataframe ##
  ##############################################
  
  step_mean <- data.frame('IS.a.Rsq'=radj_estimates[11,1], 'IS.RMSE'=rsme_estimates[11,1], 
                          'OS.a.Rsq'=radj_estimates[11,3], 'OS.RMSE'=rsme_estimates[11,3])
  boot_mean <- data.frame('IS.a.Rsq'=radj_estimates[11,2], 'IS.RMSE'=rsme_estimates[11,2], 
                          'OS.a.Rsq'=radj_estimates[11,4], 'OS.RMSE'=rsme_estimates[11,4])
  mean_est <- bind_rows(mean_est, step_mean, boot_mean)
}


############################################
## Update rownames and print final tables ##
############################################

resultrownames <- c("FC-246/0-BSR", "FC-246/0-BeEN", "FC-184/62-BSR", "FC-184/62-BeEN",
                    "IC-126/0-BSR", "IC-126/0-BeEN", "IC-94/32-BSR", "IC-94/32-BeEN")
row.names(selected_pred) <- resultrownames
selected_pred <- round(selected_pred * 100, 0)
row.names(mean_est) <- resultrownames


########################################################
## Visualize a-rsq and rsme heterogeneity in one plot ##
########################################################

grid.arrange(Full_radj, Intervention_radj, Full_rsme, Intervention_rsme, ncol=2, nrow=2)








