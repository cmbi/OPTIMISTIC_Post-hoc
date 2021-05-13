# Script run with R version 4.0.5 (2021-03-31)
# Last changes applied on: 30/04/2021

## purpose of analysis: sensitivity analysis for delta-DM1-Activ-c predictors; repeating the main analysis for the control group
## dependent variable: delta-DM1-Activ-C (10M vs baseline)
## only control group included; as such no CBT related variables considered
## excluding of patients with maximum baseline score of delta-DM1-Acitv-c
## in addition to primary/secondary outcome variables, patient/genetic variables, effect modifiers and treatment related variables are included

RNGkind(sample.kind = "Rounding") #conforms seed based sampling to earlier R versions
set.seed(1) # Allows reproducability

###############
## Libraries ##
###############

library("readxl") #V 1.3.1
library("dplyr")  #V 1.0.5
library("glmnet") #V 4.1.1
library("mice")   #V 3.13.0
library("boot")   #V 1.3.27
library("gvlma")  #V 1.0.0.3
library("car")    #V 3.0.10

###############
## Functions ##
###############

# Function to run cv.glmnet on bootstrap distributions and return coefficients
fitglmnet <- function(data, indices, alpha, nfolds){
  
  # Allow boot to select sample
  d <- data[indices,]
  
  # Matrix conversion & splitting dependent variable
  x <- as.matrix(d[,-which(names(d)=="ddm1a")])
  y <- as.matrix(d[,"ddm1a", drop=F])
  
  # Fit model and return coefficients
  fit <- cv.glmnet(x, y, type.measure ="mse", alpha = alpha, family="gaussian", nfolds = nfolds)
  coef <- coef(fit, s="lambda.1se")
  return(as.vector(coef))
}

# Function to calculate r^2 and adjusted r^2 on testingset / prediction
rsq <- function(observed, predicted, p, N){
  # p = number of predictors excluding constant; N = total sample size
  # r^2
  SSRes <- sum((observed - predicted)^2)
  SSTot <- sum((observed - mean(observed))^2)
  Rsq <- round(1 - (SSRes / SSTot),2)
  
  # Adjusted r^2
  Rsq_adj <- round(1 - (((1-Rsq) * (N - 1)) / (N - p - 1)),2)
  
  return(c(Rsq, Rsq_adj))
}

############################################
## Set working directory and open datbase ##
############################################

setwd("Z:/documents")
df <- read_excel("Outcome measures overview (V5)_corrected.xlsx", sheet = 1)


###########################
## Variable manipulation ##
###########################

## The stroop test resulted in 6 variables: stroop-time 1-3 and stroop-error 1-3
# Calculation of final stroop-interference scores is based on paper: https://doi.org/10.1016/j.nicl.2018.101615
# Accuarcy = #correct / #total; accuracy-speed trade-off time [ASTOT]: accuracy / response time 
# Stroop interference scores: ASTOT3 / ASTOT2

# stroop2 interference score at V2
stroop2V2 <- ((60 - df$StroopCardIIErrorsV2) / 60) / df$StroopCardIITimeV2
stroop3V2 <- ((60 - df$StroopCardIIIErrorsV2) / 60) / df$StroopCardIIITimeV2
StroopInterferenceV2 <- stroop3V2 / stroop2V2

## Trail-making-test: Calculation similar to the stroop test interference scores [TMTB/TMTA]
TMTV2 <- df$TMTBV2 / df$TMTAV2

## delta-dm1-activ-c scores
ddm1a <- df$DM1ActivCV4 - df$DM1ActivCV2


################################################################################
## Creating control dataframe with baseline predictors and dependent variable ##
################################################################################

# The following baselinevariables are included in the analysis:
# Primary outcome variable: DM1ActivCV2
# Secondary outcome variables: SMWTV2, PreBORGV2, AEScScoreV2, FDSSV2, MDHIV2, CISFatigueV2, BDIFsV2, INQOLQolScoreV2, stroopV2, MeanENMOV2, M5ENMOV2, L5ENMOV2
# Patient variables: AgeV2, SexCode, AgeAtOnset,  
# Genetic variables: VariantRepeats, ePAL, V2Mode, CTGDiagnostic 
# Effect modifiers: MIRSV2, TMTV2, McGIllPainV2, ASBQV2, SIPScoreV2, SSLDScoreV2, SSLNScoreV2, SSLIScoreV2, JFCSV2, ICQV2,
# IMQV2, SES28V2, CSIV2, AESIV2, CISactivityv2
# Dependent variable: ddm1a

# note. V2Diff is dropped because the variable is a linear combination of V2Mode - ePAL
bdf <- data.frame(
  df$DM1ActivCV2, df$SMWTV2, df$PreBORGV2, df$AEScScoreV2, df$FDSSV2, df$MDHIV2, df$CISFatigueV2, df$BDIFsV2, df$INQOLQolScoreV2, StroopInterferenceV2, df$MeanENMOV2, df$M5ENMOV2, df$L5ENMOV2,
  df$AgeV2, df$SexCode, df$AgeAtOnset, df$VariantRepeats,
  df$ePAL, df$V2Mode, df$CTGDiagnostic,
  df$MIRSV2, TMTV2, df$McGillPainV2, df$ASBQV2, df$SIPScoreV2, df$SSLDScoreV2, df$SSLNScoreV2, df$SSLIScoreV2, df$JFCSV2, df$ICQV2,
  df$IMQV2, df$SES28V2, df$CSIV2, df$AESIV2, df$CISactivityV2,
  ddm1a)

# Clean up variable names & change all variable types to numeric
colnames(bdf) <- gsub(x = colnames(bdf), pattern = "df.", replacement = "")
colnames(bdf) <- gsub(x = colnames(bdf), pattern = "V2", replacement = "")
bdf <- sapply(bdf[,1:length(colnames(bdf))], function(x){as.numeric(as.character(x))})
bdf <- as.data.frame(bdf)

# Subset the control group #n=127
bdf <- bdf[df$TreatmentCode==0,] 


#######################
## Handling outliers ##
#######################

## Screening for outliers
# Outliers defined as < Q1-2*IQR | > Q3+2*IQR

cdf <- bdf #a copy is made in order to detect differnces after removal
for (x in 1:length(colnames(bdf))){
  Q1 <- unname(quantile(bdf[,x], 0.25, na.rm = T))
  Q3 <- unname(quantile(bdf[,x], 0.75, na.rm = T))
  IQR <- Q3-Q1
  l <- Q1 - 2*IQR
  u <- Q3 + 2*IQR
  bdf[,x][bdf[,x] > u | bdf[,x] < l] <- NA
}

# Number of values changed per column (=differences in NA values)
a <- apply(cdf, 2, function(x) length(which(!is.na(x)))) - apply(bdf, 2, function(x) length(which(!is.na(x))))
sum(apply(cdf, 2, function(x) length(which(!is.na(x)))) - apply(bdf, 2, function(x) length(which(!is.na(x)))))


cdf$StroopInterference[cdf$StroopInterference > 1.1] <- NA #1
cdf$M5ENMO[cdf$M5ENMO > 200] <- NA #1
cdf$L5ENMO[cdf$L5ENMO > 17] <- NA #2
cdf$SSLNScore[cdf$SSLNScore > 20] <- NA #1
cdf$ddm1a[cdf$ddm1a < -70] <- NA #1


# restore dataframe - 6 outliers removed
bdf <- cdf
cdf <- NULL


############################
## Missing Value analysis ##
############################

# Exclude all patients that have no ddm1a score (can neither be used for training nor testing)
bdf <- bdf[!is.na(bdf$ddm1a),] # Note: 115 cases remain (12 cases removed)

# Exclude all patients that have max baseline DM1-Activ-C scores [7 patients removed - 108 remain]
bdf <- bdf[!bdf$DM1ActivC==100,]

# Number of NA's per column
na_count <- sapply(bdf, function(y) sum(length(which(is.na(y)))))
na_count <- data.frame(na_count)
na_count["percentage"] <- round((na_count/nrow(bdf))*100, 2)
na_count[order(-na_count$percentage),]


# variables with > 25% missing values are removed
# Excluding SIPScore (n=77, 71,30%), AESI (n=64, 59,26%), CSI (n = 63, 58,33%), CTGDiagnostic (n=57, 52,78%)

bdf$SIPScore <- NULL
bdf$CSI <- NULL
bdf$AESI <- NULL
bdf$CTGDiagnostic <- NULL


##################################
## Imputation of missing values ##
##################################

imputed <- mice(bdf, m=10, maxit=50, meth= 'pmm', print=T) 


#############################################
## Predictor selection on imputed datasets ##
#############################################

## Dataframe to store selected predictors
selected_pred <- data.frame(matrix(ncol=32, nrow=0))
colnames(selected_pred) <- c("Intercept", colnames(bdf)[-32])

for (i in 1:10){
  ## Select imputed dataset
  data <- complete(imputed, i)
  
  ## Fit Elastic-Net on 5000 bootstrap distributions | 5 fold CV, alpha=0.5
  boot_res <- boot(data=data, fitglmnet, R=5000, alpha=0.5, nfolds=5)
  boot_freq <- as.data.frame(boot_res$t)
  colnames(boot_freq) <- colnames(selected_pred)
  
  ## Calculate variable inclusion frequencies
  boot_sums <- round((colSums(boot_freq != 0) / 5000), 2)
  
  ## Store results
  selected_pred <- bind_rows(selected_pred, boot_sums)
}


################
## Store VIPs ##
################

## Select imputation independent frequently selected variables | inclusion frequency >= 60%
cm <- colMeans(selected_pred)
predictors <- names(cm[cm >= 0.6])
predictors <- predictors[predictors != "Intercept"]

## Add mean frequencies to predictor dataframe & print+store ordered table
selected_pred <- rbind(selected_pred, cm)
row.names(selected_pred) <- c(1:10, "Mean")
t_selected_pred <- as.data.frame(t(selected_pred))
row.names(t_selected_pred) <- gsub(x = row.names(t_selected_pred), pattern = "SMWT", replacement = "6MWT")
t_selected_pred[order(-t_selected_pred$Mean),]
t_selected_pred <- t_selected_pred * 100
t_selected_pred$Mean <- round(t_selected_pred$Mean, 0)

write.csv(t_selected_pred[order(-t_selected_pred$Mean),], "OPH_dDM1Activc_senscontrol_VIPs.csv", row.names = TRUE)
write.csv(t_selected_pred[order(-t_selected_pred$Mean),], "OPH_dDM1Activc_senscontrol_VIPs.xlsx", row.names = TRUE)





