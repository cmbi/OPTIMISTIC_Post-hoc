# Script run with R version 4.0.5 (2021-03-31)
# Last changes applied on: 30/04/2021

## purpose of analysis: identification of variables associated with CBT treatment
## dependent variable: delta-6MWT (10M vs baseline)
## only treatment group included
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
  x <- as.matrix(d[,-which(names(d)=="d6MWT")])
  y <- as.matrix(d[,"d6MWT", drop=F])
  
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
df2 <- read_excel("Outcome measures overview (V5)_corrected.xlsx", sheet = 3)


######################
## DF2 manipulation ##
######################

## Subsetting dataframe
df2 <- data.frame(df2$`Randomisation number`, df2$`Duration of session`, df2$`Session number`,
                  df2$`Module dealt with in session: module 1 (goal-setting)`,df2$`Module dealt with in session: module 2 (sleep-wake pattern)`, df2$`Module dealt with in session: module 3 (getting started)`,
                  df2$`Module dealt with in session: module 4 (Activity)`, df2$`Module dealt with in session: module 5 (Helpful beliefs)`, df2$`Module dealt with in session: module 6 (Optimizing social interactions)`,
                  df2$`Module dealt with in session: module 7 (Pain)`,
                  df2$`Modules indicated as per baseline assessment? Module 1`, df2$`Modules indicated as per baseline assessment? Module 2`, df2$`Modules indicated as per baseline assessment? Module 3`,
                  df2$`Modules indicated as per baseline assessment? Module 4`, df2$`Modules indicated as per baseline assessment? Module 5`, df2$`Modules indicated as per baseline assessment? Module 6`,
                  df2$`Modules indicated as per baseline assessment? Module 7`)

# Change relevant variables to type 'numeric'
df3 <- as.data.frame(sapply(df2[,2:length(colnames(df2))], function(x){as.numeric(as.character(x))}))
df3['Randomisaton number'] <- as.character(df2$df2..Randomisation.number.)
colnames(df3) <- gsub(x = colnames(df3), pattern = "df2..", replacement = "")
df2 <- df3
df3 <- NULL

## Imputation of session time: average of the other sessions
df2[609, 1] <- mean(df2[df2$`Randomisaton number`=='NCL056P',1], na.rm=T)
A008P <- mean(df2[df2$`Randomisaton number`=='A008P',1], na.rm=T)
df2[675, 1] <- A008P
df2[676, 1] <- A008P

## Dropping 2 patients as they didnt recieve therapy
df2 <- df2[!df2$`Randomisaton number`=='M053P',]
df2 <- df2[!df2$`Randomisaton number`=='M061P',]

## Dropping 2 patients because of > 1 indications missing
df2 <- df2[!df2$`Randomisaton number`=='NCL002P',]
df2 <- df2[!df2$`Randomisaton number`=='A078P',]


##########################
## variable calculation ##
##########################

## indication score [% of how many indicated modules were atleast given once]

# mdf stores which modules were at least once given
mdf <- data.frame(df2$Module.dealt.with.in.session..module.1..goal.setting.., df2$Module.dealt.with.in.session..module.2..sleep.wake.pattern.., df2$Module.dealt.with.in.session..module.3..getting.started..,
                  df2$Module.dealt.with.in.session..module.4..Activity..,df2$Module.dealt.with.in.session..module.5..Helpful.beliefs..,df2$Module.dealt.with.in.session..module.6..Optimizing.social.interactions..,
                  df2$Module.dealt.with.in.session..module.7..Pain..)
mdf[is.na(mdf)] <- 0 #optional, however no missing values expected
mdf <- aggregate(x = mdf,
                 by = list(RN = df2$`Randomisaton number`),
                 FUN = sum)
RN <- mdf$RN
mdf[mdf > 1] <- 1
mdf$RN <- RN

# idf stores which modules were indicated
idf <- data.frame(df2$Modules.indicated.as.per.baseline.assessment..Module.1., df2$Modules.indicated.as.per.baseline.assessment..Module.2.,df2$Modules.indicated.as.per.baseline.assessment..Module.3.,
                  df2$Modules.indicated.as.per.baseline.assessment..Module.4., df2$Modules.indicated.as.per.baseline.assessment..Module.5., df2$Modules.indicated.as.per.baseline.assessment..Module.6.,
                  df2$Modules.indicated.as.per.baseline.assessment..Module.7.)
idf[is.na(idf)] <- 0 
idf <- aggregate(x = idf,
               by = list(RN = df2$`Randomisaton number`),
               FUN = sum)

# store total number of indications per patient
idf$sum <- rowSums(idf[,2:length(idf)])

# Calculate per patient and module whether or not is has been completed at least once
m1 <- mdf$df2.Module.dealt.with.in.session..module.1..goal.setting../idf$df2.Modules.indicated.as.per.baseline.assessment..Module.1.
m2 <- mdf$df2.Module.dealt.with.in.session..module.2..sleep.wake.pattern../idf$df2.Modules.indicated.as.per.baseline.assessment..Module.2.
m3 <- mdf$df2.Module.dealt.with.in.session..module.3..getting.started../idf$df2.Modules.indicated.as.per.baseline.assessment..Module.3.
m4 <- mdf$df2.Module.dealt.with.in.session..module.4..Activity../idf$df2.Modules.indicated.as.per.baseline.assessment..Module.4.
m5 <- mdf$df2.Module.dealt.with.in.session..module.5..Helpful.beliefs../idf$df2.Modules.indicated.as.per.baseline.assessment..Module.5.
m6 <- mdf$df2.Module.dealt.with.in.session..module.6..Optimizing.social.interactions../idf$df2.Modules.indicated.as.per.baseline.assessment..Module.6.
m7 <- mdf$df2.Module.dealt.with.in.session..module.7..Pain../idf$df2.Modules.indicated.as.per.baseline.assessment..Module.7.

complete <- data.frame(m1, m2,m3,m4,m5,m6,m7)
complete[complete == Inf] <- 0 # inf means: module completed that was not indicated
complete[is.na(complete)] <- 0 # NaN means: not completed module that was not indicated 
complete <- data.frame(RN=mdf$RN, nIndications=idf$sum, IndScore=round(rowSums(complete)/idf$sum,4))


## Number of sessions & Total session duration
tdf <- data.frame(nSessions=df2$Session.number., Sessiontime=df2$Duration.of.session.)
tdf$nSessions <- 1 #sessions are numbered, by assigning each "1" it is possible to count sessions easily
tdf <- aggregate(x = tdf,
                 by = list(RN = df2$`Randomisaton number`),
                 FUN = sum)

## Storage in completed dataframe
cdf2 <- merge(tdf, complete, by='RN')

## Final check for completeness
cdf2 <- cdf2[!cdf2$RN =='M066P',] # missing session time and modules dealt with; only 1 session completed
## complete treatment information dataframe n=121


############################
## Combinging df1 and df2 ##
############################

## Subgrouping intervention arm
df <- data.frame(df)
df <- df[df$TreatmentCode=='1',]

## Nomenclature missmatch correction
cdf2$RN <- gsub(x=cdf2$RN, pattern="NCL", replacement ="D")
cdf2$RN <- gsub(x=cdf2$RN, pattern="M", replacement ="B")

# Removing patients in df that have been dropped in df2: M053P (B053P), M061P (B061P), NCL002P (D002P), M066P(B066P)
df <- df[!df$PatientID...1=='B053P',]
df <- df[!df$PatientID...1=='B061P',]
df <- df[!df$PatientID...1=='D002P',]
df <- df[!df$PatientID...1=='A078P',]
df <- df[!df$PatientID...1=='B066P',]

# Comparing patient ID's
names1 <- df$PatientID...1
names2 <- cdf2$RN
setdiff(names1, names2)

# Dropping two more patients in df that are not present in cdf2
df <- df[!df$PatientID...1=='A027P',]
df <- df[!df$PatientID...1=='B025P',]

# Merging the two dataframes based on patientID
colnames(df)[1] <- "RN"
df <- merge(df, cdf2, by='RN')


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

## delta-6MWT scores
d6MWT <- df$SMWTV4-df$SMWTV2


########################################################################
## Creating dataframe with baseline predictors and dependent variable ##
########################################################################

# The following baselinevariables are included in the analysis:
# Primary outcome variable: DM1ActivCV2
# Secondary outcome variables: SMWTV2, PreBORGV2, AEScScoreV2, FDSSV2, MDHIV2, CISFatigueV2, BDIFsV2, INQOLQolScoreV2, stroopV2, MeanENMOV2, M5ENMOV2, L5ENMOV2
# Patient variables: AgeV2, SexCode, AgeAtOnset,  
# Genetic variables: VariantRepeats, ePAL, V2Mode, CTGDiagnostic 
# Effect modifiers: GradedExerciseTherapy, MIRSV2, TMTV2, McGIllPainV2, ASBQV2, SIPScoreV2, SSLDScoreV2, SSLNScoreV2, SSLIScoreV2, JFCSV2, ICQV2,
# IMQV2, SES28V2, CSIV2, AESIV2, CISactivityv2
# Therapy variables: nsessions, sessiontime, indications, indscore
# Dependent variable: d6MWT

# note. V2Diff is dropped because the variable is a linear combination of V2Mode - ePAL
bdf <- data.frame(
  df$DM1ActivCV2, df$SMWTV2, df$PreBORGV2, df$AEScScoreV2, df$FDSSV2, df$MDHIV2, df$CISFatigueV2, df$BDIFsV2, df$INQOLQolScoreV2, StroopInterferenceV2, df$MeanENMOV2, df$M5ENMOV2, df$L5ENMOV2,
  df$AgeV2, df$SexCode, df$AgeAtOnset, df$VariantRepeats,
  df$ePAL, df$V2Mode, df$CTGDiagnostic,
  df$GradedExerciseTherapy, df$MIRSV2, TMTV2, df$McGillPainV2, df$ASBQV2, df$SIPScoreV2, df$SSLDScoreV2, df$SSLNScoreV2, df$SSLIScoreV2, df$JFCSV2, df$ICQV2,
  df$IMQV2, df$SES28V2, df$CSIV2, df$AESIV2, df$CISactivityV2,
  df$nSessions, df$Sessiontime, df$nIndications, df$IndScore,
  d6MWT)

# Clean up variable names & change all variable types to numeric
colnames(bdf) <- gsub(x = colnames(bdf), pattern = "df.", replacement = "")
colnames(bdf) <- gsub(x = colnames(bdf), pattern = "V2", replacement = "")
bdf <- sapply(bdf[,1:length(colnames(bdf))], function(x){as.numeric(as.character(x))})
bdf <- as.data.frame(bdf)


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

cdf$TMT[cdf$TMT > 5] <- NA #1
cdf$SSLNScore[cdf$SSLNScore > 25] <- NA #1
cdf$d6MWT[cdf$d6MWT < -200 | cdf$d6MWT > 200] <- NA#2

# restore dataframe - 4 outliers removed
bdf <- cdf
cdf <- NULL


############################
## Missing Value analysis ##
############################

# Exclude all patients that have no d6MWT score (can neither be used for training nor testing)
bdf <- bdf[!is.na(bdf$d6MWT),] # Note: 103 cases remain (18 cases removed)

# Number of NA's per column
na_count <- sapply(bdf, function(y) sum(length(which(is.na(y)))))
na_count <- data.frame(na_count)
na_count["percentage"] <- round((na_count/nrow(bdf))*100, 2)

row.names(na_count) <- gsub(x = row.names(na_count), pattern = "SMWT", replacement = "6MWT")
na_count[order(-na_count$percentage),]


# variables with > 25% missing values are removed
# Excluding CSI(n=61, 59,22%), AESI (n=61, 59,22%), CTGDiagnostic (n=58, 56,31%)
bdf$CTGDiagnostic <- NULL
bdf$CSI <- NULL
bdf$AESI <- NULL


##################################
## Imputation of missing values ##
##################################

imputed <- mice(bdf, m=10, maxit=50, meth= 'pmm', print=T) 
#densityplot(imputed)
#plot(imputed, c("MeanENMO", "M5ENMO", "L5ENMO"))
#plot(imputed, c("ASBQ", "SSLDScore", "AgeAtOnset"))
#plot(imputed, c("SSLIScore", "McGillPain"))


#############################################
## Predictor selection on imputed datasets ##
#############################################

## Dataframe to store selected predictors
selected_pred <- data.frame(matrix(ncol=38, nrow=0))
colnames(selected_pred) <- c("Intercept", colnames(bdf)[-38])

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


########################
## Build linear model ##
########################

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

write.csv(t_selected_pred[order(-t_selected_pred$Mean),], "OPH_6MWT_VIPs.csv", row.names = TRUE)
write.csv(t_selected_pred[order(-t_selected_pred$Mean),], "OPH_6MWT.xlsx", row.names = TRUE)

## Create regression formula
formula <- paste(predictors, collapse = " + ")
formula <- paste("d6MWT", formula, sep=" ~ ")

## Build model & check linear regression assumptions
fit <- lm(as.formula(formula), data=bdf)
gvlma(fit)
summary(fit)

## Detect and remove signifcant outliers
outliers <- outlierTest(fit)
outliers <- c("54")
bdf <- bdf[!(row.names(bdf) %in% outliers),]

## Build model & check lienar regression assumptions
fit <- lm(as.formula(formula), data=bdf)
summary(fit)
gvlma(fit)
rmse <- round(sqrt(mean(fit$residuals^2)),2)



