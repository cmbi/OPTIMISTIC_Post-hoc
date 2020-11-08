# Script run with R version 3.5.1 (2018-07-12)
# Last changes applied on: 30/07/2020

# Script to calculate spearman rank's correlation coefficients among research measurements used in the OPTIMISTIC clinical trial
# Significance assessed by pairwise t-tests based on spearman rank's rho, FDR set to 5% using Benjamini Hochberg


###############
## Libraries ##
###############
library("readxl")     #1.3.1
library("corrplot")   #0.84
library("psych")      #1.8.12
library("ggplot2")    #3.1.1
library("gridExtra")  #2.3


############################################
## Set working directory and open databse ##
############################################
setwd("Z:/documents")
df <- read_excel("Outcome measures overview (V5)_corrected.xlsx", sheet = 1)


###########################
## Variable manipulation ##
###########################

## The stroop test resulted in 6 variables: stroop-time 1-3 and stroop-error 1-3
# Calculation based on publications, see table 1 in paper for references
# accuarcy = #correct / #total; accuracy-speed trade-off time [ASTOT]: accuracy / response time 
# stroop interference scores: ASTOT3 / ASTOT2

# stroop2 interference score at V2
stroop2V2 <- ((60 - df$StroopCardIIErrorsV2) / 60) / df$StroopCardIITimeV2
stroop3V2 <- ((60 - df$StroopCardIIIErrorsV2) / 60) / df$StroopCardIIITimeV2
StroopInterferenceV2 <- stroop3V2 / stroop2V2

# stroop interference score at V3
stroop2V3 <- ((60 - df$StroopCardIIErrorsV3) / 60) / df$StroopCardIITimeV3
stroop3V3 <- ((60 - df$StroopCardIIIErrorsV3) / 60) / df$StroopCardIIITimeV3
dStroopInterferenceV3 <- stroop3V3 - stroop2V3

# Stroop interference scores at V4
stroop2V4 <- ((60 - df$StroopCardIIErrorsV4) / 60) / df$StroopCardIITimeV4
stroop3V4 <- ((60 - df$StroopCardIIIErrorsV4) / 60) / df$StroopCardIIITimeV4
dStroopInterferenceV4 <- stroop3V4 / stroop2V4

# Stroop interference scores at V5
stroop2V5 <- ((60 - df$StroopCardIIErrorsV5) / 60) / df$StroopCardIITimeV5
stroop3V5 <- ((60 - df$StroopCardIIIErrorsV5) / 60) / df$StroopCardIIITimeV5
dStroopInterferenceV5 <- stroop3V5 / stroop2V5


## Trail-making-test: Calculation similar to the stroop test interference scores [TMTB/TMTA]
TMTV2 <- df$TMTBV2 / df$TMTAV2
dTMTV3 <- df$TMTBV3 / df$TMTAV3
dTMTV4 <- df$TMTBV4 / df$TMTAV4
dTMTV5 <- df$TMTBV5 / df$TMTAV5


############################
## Correlation analysis 1 ##
############################
# Correlation of baseline (t=0 measurements)

## Subset all relevant baseline outcome measurements
# The following baselinevariables are included in the analysis:
# Primary outcome variable: DM1ActivCV2
# Secondary outcome variables: 6MWTV2, PreBORGV2, AEScScoreV2, FDSSV2, MDHIV2, CISFatigueV2, BDIFsV2, INQOLQolScoreV2, stroopV2, MeanENMOV2, M5ENMOV2, L5ENMOV2
# Effect modifiers: MIRSV2, TMTV2, McGIllPainV2, ASBQV2, SSLDScoreV2, SSLNScoreV2, SSLIScoreV2, JFCSV2, ICQV2,
# IMQV2, SES28V2, CSIV2, AESIV2, CISactivityv2

INQoLV2 <- df$INQOLQolScoreV2
BDIFSV2 <- df$BDIFsV2

bdf <- data.frame(
  df$DM1ActivCV2, 
  df$SMWTV2, df$PreBORGV2, df$AEScScoreV2, df$FDSSV2, df$MDHIV2, df$CISFatigueV2, BDIFSV2, INQoLV2, StroopInterferenceV2, df$MeanENMOV2, df$M5ENMOV2, df$L5ENMOV2,
  df$MIRSV2, TMTV2, df$McGillPainV2, df$ASBQV2, df$SSLDScoreV2, df$SSLNScoreV2, df$SSLIScoreV2, df$JFCSV2, df$ICQV2,
  df$IMQV2, df$SES28V2, df$CSIV2, df$AESIV2, df$CISactivityV2)

## Check datatype of each variable
sapply(bdf, class) # -> each variable is of type 'numeric', no conversion necessary

## Clean up variable names
colnames(bdf) <- gsub(x = colnames(bdf), pattern = "df.", replacement = "")
colnames(bdf) <- gsub(x = colnames(bdf), pattern = "V2", replacement = "")
colnames(bdf) <- gsub(x = colnames(bdf), pattern = "DM1ActivC", replacement = "DM1-Activ-c")
colnames(bdf) <- gsub(x = colnames(bdf), pattern = "SMWT", replacement = "6MWT")


## Shapiro-Wilks test for normality
pvalues <- apply(bdf, 2, function(x) shapiro.test(x)$p.value)
round(pvalues[pvalues > 0.05], 2) # 4/27 normally distributed

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
apply(cdf, 2, function(x) length(which(!is.na(x)))) - apply(bdf, 2, function(x) length(which(!is.na(x))))
sum(apply(cdf, 2, function(x) length(which(!is.na(x)))) - apply(bdf, 2, function(x) length(which(!is.na(x)))))

# 20 outliers detected, visual inspection confirmed 6 likely true outliers
cdf$StroopInterference[cdf$StroopInterference > 1.1] <- NA #1
cdf$MeanENMO[cdf$MeanENMO > 60] <- NA #2
cdf$M5ENMO[cdf$M5ENMO > 200] <- NA #1
cdf$L5ENMO[cdf$L5ENMO > 17] <- NA #2

## restoring dataset and removing 6 outliers
bdf <- cdf

## Correlation analysis
# psych package using corr.test function
# non-parametric spearman correlations calculated given the non-normally distributed data
# benjamini hochberg correction applied for multiple testing
# for symmetric matrices, raw probabilites are reported below the diagonal and correlations adjusted for multiple comparisons above the diagonal

bcordfbh <- corr.test(bdf, use = "pairwise",method="spearman",adjust="BH", 
                      alpha=.05)


############################
## Correlation analysis 2 ##
############################

## Subset all relevant delta outcome measurements at t=10M
# The following primary endpoint variables are included in the analysis:
# Primary outcome variable: DM1ActivCV4
# Secondary outcome variables: 6MWTV4, PreBORGV4, AEScScoreV4, FDSSV4, MDHIV4, CISFatigueV4, BDIFsV4, INQOLQolScoreV4, stroopV4, MeanENMOV4, M5ENMOV4, L5ENMOV4
# Effect modifiers: MIRSV4, TMTV4, McGIllPainV4, ASBQV4, SSLDScoreV4, SSLNScoreV4, SSLIScoreV4, JFCSV4, ICQV4,
# IMQV4, SES28V4, CSIV4, AESIV4, CISactivityV4

# These variables are added seperately to adjust their names individually
dINQoLV4 <- df$INQOLQolScoreV4
dBDIFSV4 <- df$BDIFsV4
dMeanENMOV4 <- df$`MeanENMO  V4`
dM5ENMOV4 <- df$`M5ENMO  V4`

V4 <- data.frame(
  df$DM1ActivCV4, 
  df$SMWTV4, df$PreBORGV4, df$AEScScoreV4, df$FDSSV4, df$MDHIV4, df$CISFatigueV4, dBDIFSV4, dINQoLV4, dStroopInterferenceV4, dMeanENMOV4, dM5ENMOV4, df$L5ENMOV4,
  df$MIRSV4, dTMTV4, df$McGillPainV4, df$ASBQV4, df$SSLDScoreV4, df$SSLNScoreV4, df$SSLIScoreV4, df$JFCSV4, df$ICQV4,
  df$IMQV4, df$SES28V4, df$CSIV4, df$AESIV4, df$CISactivityV4)

# Check datatype of each variable
sapply(V4, class) # -> each variable is of type 'numeric', no conversion necessary

## Calculate delta-dataframe by subtracting bdf from edf
ddf <- V4 - bdf

## Clean up variable names
colnames(ddf) <- gsub(x = colnames(ddf), pattern = "df.", replacement = "d")
colnames(ddf) <- gsub(x = colnames(ddf), pattern = "V4", replacement = "")
colnames(ddf) <- gsub(x = colnames(ddf), pattern = "dDM1ActivC", replacement = "dDM1-Activ-c")
colnames(ddf) <- gsub(x = colnames(ddf), pattern = "dSMWT", replacement = "d6MWT")

## Restrict corrleation analysis to intervention group
ddf <- ddf[df$Treatment == "Intervention",] #n=128

## Shapiro-Wilks test for normality
pvalues <- apply(ddf, 2, function(x) shapiro.test(x)$p.value)
round(pvalues[pvalues > 0.05], 2) # 8/27 normally distributed

## Screening for outliers
# Outliers defined as < Q1-2*IQR | > Q3+2*IQR

cdf <- ddf #a copy is made in order to detect differnces after removal
for (x in 1:length(colnames(ddf))){
  Q1 <- unname(quantile(ddf[,x], 0.25, na.rm = T))
  Q3 <- unname(quantile(ddf[,x], 0.75, na.rm = T))
  IQR <- Q3-Q1
  l <- Q1 - 2*IQR
  u <- Q3 + 2*IQR
  ddf[,x][ddf[,x] > u | ddf[,x] < l] <- NA
}

# Number of values changed per column (=differences in NA values)
apply(cdf, 2, function(x) length(which(!is.na(x)))) - apply(ddf, 2, function(x) length(which(!is.na(x))))
sum(apply(cdf, 2, function(x) length(which(!is.na(x)))) - apply(ddf, 2, function(x) length(which(!is.na(x)))))
# -> 96 potential outliers identified, 18 removed

cdf$d6MWT[cdf$d6MWT < -200 | cdf$d6MWT > 300] <- NA #2
cdf$dMDHI[cdf$dMDHI < -40 | cdf$dMDHI > 30] <- NA #3
cdf$dBDIFS[cdf$dBDIFS > 9] <-  NA #1
cdf$dTMT[cdf$dTMT > 4] <- NA #1
cdf$dMcGillPain[cdf$dMcGillPain > 50] <- NA #1
cdf$dASBQ[cdf$dASBQ > 16] <- NA #3
cdf$dSSLNScore[cdf$dSSLNScore > 10] <- NA #1
cdf$dSSLIScore[cdf$dSSLIScore < -40] <- NA #1
cdf$dICQ[cdf$dICQ < -10] <- NA #1
cdf$dIMQ[cdf$dIMQ < -30] <- NA #1
cdf$dSES28[cdf$dSES28 > 10] <- NA #2
cdf$dCSI[cdf$dCSI > 5] <- NA #1

## restoring dataset and removing 18 outliers
ddf <- cdf

## Correlation analysis
# psych package using corr.test function
# non-parametric spearman correlations calculated given the non-normally distributed data
# Benjamini Hochberg correction applied for multiple testing
# for symmetric matrices, raw probabilites are reported below the diagonal and correlations adjusted for multiple comparisons above the diagonal
dcordfbh <- corr.test(ddf, use ="pairwise", method="spearman", adjust="BH", alpha = .05)


########################################
## Visualization of corr 1 and corr 2 ##
########################################

plot.new()
tiff("OPH_Fig1_Correlograms.tiff", units="cm", width=13, height=21.7, res=300)
par(mfrow=c(2,1))
corrplot(bcordfbh$r, tl.col = "black", tl.srt = 45, tl.cex = 0.60, cl.cex = 0.60, mar=c(0,0,0,0),
         p.mat = bcordfbh$p, sig.level = 0.05, insig = "blank", family="sans")
corrplot(dcordfbh$r, tl.col = "black", tl.srt = 45, tl.cex = 0.60, cl.cex = 0.60, mar=c(0,0,0,0),
         p.mat = dcordfbh$p, sig.level = 0.05, insig = "blank", family="sans")
mtext("A", side = 1, line = -36.5, cex=1, adj=0.07, font=2)
mtext("B", side = 1, line = -15.5, cex=1, adj=0.07, font=2)
dev.off()














