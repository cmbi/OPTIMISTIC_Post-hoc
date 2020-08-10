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

par(mfrow=c(1,2))
corrplot(bcordfbh$r, tl.col = "black", tl.srt = 45, 
         p.mat = bcordfbh$p, sig.level = 0.05, insig = "blank")
mtext("A", side = 1, line = -45, cex=2, adj=0.05, font=2)

corrplot(dcordfbh$r, tl.col = "black", tl.srt = 45, 
         p.mat = dcordfbh$p, sig.level = 0.05, insig = "blank")
mtext("B", side = 1, line = -45, cex=2, adj=0.09, font=2)

# Exported using RStudio as 1900*1200 JPEG

##########################################
## Scatter plots of DM1-Activ-c vs 6MWT ##
##########################################

## Scatterplot of baseline DM1ActivC and 6MWT
baseplot <- ggplot(bdf, aes_string(x=bdf$`DM1-Activ-c`, y=bdf$'6MWT')) +
  geom_point(color=("black"), size=2) +
  labs(x="DM1-Activ-c (t = 0M)", y="6MWT (t = 0M)", tag="A") +
  scale_x_continuous(expand = c(0,0), limits = c(20, 105), breaks = scales::pretty_breaks(n=10)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 800), breaks = scales::pretty_breaks(n=10)) +
  geom_smooth(method="lm", color=("black"), se=F) +
  theme(
    panel.border = element_rect(colour="black", fill = NA, size = 0.4),
    axis.text = element_text(size = 15),
    axis.title = element_text(colour="black", size = 20, face="bold"),
    plot.tag.position = "topleft",
    plot.tag = element_text(size=30, face="bold"),
    panel.background = element_rect(fill="white"),
    aspect.ratio = 1
  )

## Scatterplot of dDM1ActivC and d6MWT
deltaplot <- ggplot(ddf, aes_string(x=ddf$`dDM1-Activ-c`, y=ddf$d6MWT)) +
  geom_point(color=("black"), size=2) +
  labs(x="dDM1-Activ-c", y="d6MWT", tag="B") +
  scale_x_continuous(expand = c(0,0), limits = c(-25, 36), breaks = scales::pretty_breaks(n=10)) +
  scale_y_continuous(expand = c(0,0), limits = c(-125, 205), breaks = scales::pretty_breaks(n=10)) +
  geom_smooth(method="lm", color=("black"), se=F) +
  theme(
    panel.border = element_rect(colour="black", fill = NA, size = 0.4),
    axis.text = element_text(size = 15),
    axis.title = element_text(colour="black", size = 20, face="bold"),
    plot.tag.position = "topleft",
    plot.tag = element_text(size=30, face="bold"),
    panel.background = element_rect(fill="white"),
    aspect.ratio = 1
  )
grid.arrange(baseplot, deltaplot, ncol=2)

# Exported using RStudio as 1900*1200 JPEG

##############
## Appendix ##
##############

## As a measure of robustness, delta 5 months and delta 16 months correlations have also been calculated
## These results are reportes in de Appendix of the paper

###########################
## 5M-Delta-Correlations ##
###########################

## Subset all relevant baseline outcome measurements at t=5M
# The following primary endpoint variables are included in the analysis:
# Primary outcome variable: DM1ActivCV3
# Secondary outcome variables: 6MWTV3, PreBORGV3, AEScScoreV3, FDSSV3, MDHIV3, CISFatigueV3, BDIFsV3, INQOLQolScoreV3, stroopV3, MeanENMOV3, M5ENMOV3, L5ENMOV3
# Effect modifiers: MIRSV3, TMTV3, McGIllPainV3, ASBQV3, SSLDScoreV3, SSLNScoreV3, SSLIScoreV3, JFCSV3, ICQV3,
# IMQV3, SES28V3, CSIV3, AESIV3, CISactivityV3

# These variables are added seperately to adjust their names individually
df.INQoLV3 <- df$INQOLQolScoreV3
df.StroopInterferenceV3 <- dStroopInterferenceV3
df.TMTV3 <- dTMTV3

V3df <- data.frame(
  df$DM1ActivCV3, 
  df$SMWTV3, df$PreBORGV3, df$AEScScoreV3, df$FDSSV3, df$MDHIV3, df$CISFatigueV3, df$BDIFsV3, df.INQoLV3, df.StroopInterferenceV3, df$MeanENMOV3, df$M5ENMOV3, df$L5ENMOV3,
  df$MIRSV3, df.TMTV3, df$McGillPainV3, df$ASBQV3, df$SSLDScoreV3, df$SSLNScoreV3, df$SSLIScoreV3, df$JFCSV3, df$ICQV3,
  df$IMQV3, df$SES28V3, df$CSIV3, df$AESIV3, df$CISactivityV3)

# Check datatype of each variable
sapply(V3df, class) # -> each variable is of type 'numeric', no conversion necessary

## Calculate delta-dataframe by subtracting bdf from edf
d5df <- V3df - bdf

## Clean up variable names
colnames(d5df) <- gsub(x = colnames(d5df), pattern = "df.", replacement = "d5")
colnames(d5df) <- gsub(x = colnames(d5df), pattern = "V3", replacement = "")
colnames(d5df) <- gsub(x = colnames(d5df), pattern = "d5DM1ActivC", replacement = "d5DM1-Activ-c")
colnames(d5df) <- gsub(x = colnames(d5df), pattern = "d5SMWT", replacement = "d56MWT")

## Restrict corrleation analysis to intervention group
d5df <- d5df[df$Treatment == "Intervention",] #n=128

## Shapiro-Wilks test for normality
pvalues <- apply(d5df, 2, function(x) shapiro.test(x)$p.value)
round(pvalues[pvalues > 0.05], 2) # 9/27 normally distributed

## Screening for outliers
# Outliers defined as < Q1-2*IQR | > Q3+2*IQR

cdf <- d5df #a copy is made in order to detect differnces after removal
for (x in 1:length(colnames(d5df))){
  Q1 <- unname(quantile(d5df[,x], 0.25, na.rm = T))
  Q3 <- unname(quantile(d5df[,x], 0.75, na.rm = T))
  IQR <- Q3-Q1
  l <- Q1 - 2*IQR
  u <- Q3 + 2*IQR
  d5df[,x][d5df[,x] > u | d5df[,x] < l] <- NA
}

# Number of values changed per column (=differences in NA values)
apply(cdf, 2, function(x) length(which(!is.na(x)))) - apply(d5df, 2, function(x) length(which(!is.na(x))))
sum(apply(cdf, 2, function(x) length(which(!is.na(x)))) - apply(d5df, 2, function(x) length(which(!is.na(x)))))
# -> 89 potential outliers identified, 17 removed

cdf$`d5DM1-Activ-c`[cdf$`d5DM1-Activ-c` < -55] <- NA #1
cdf$d56MWT[cdf$d56MWT > 150] <- NA #2
cdf$d5MDHI[cdf$d5MDHI < -50 | cdf$d5MDHI > 30] <- NA #4
cdf$d5CISFatigue[cdf$d5CISFatigue > 27] <- NA #1
cdf$d5MeanENMO[cdf$d5MeanENMO < -20 | cdf$d5MeanENMO > 15] <- NA #3
cdf$d5L5ENMO[cdf$d5L5ENMO > 3] <- NA #3
cdf$d5SSLNScore[cdf$d5SSLNScore > 14] <- NA #1
cdf$d5SSLIScore[cdf$d5SSLIScore < -40] <- NA #1
cdf$d5JFCS[cdf$d5JFCS > 20] <- NA #1
cdf$d5IMQ[cdf$d5IMQ > 25] <- NA #1
cdf$d5AESI[cdf$d5AESI < -189] <- NA #1

## restoring dataset and removing 41 outliers
d5df <- cdf

## Correlation analysis
# psych package using corr.test function
# non-parametric spearman correlations calculated given the non-normally distributed data
# Benjamini Hochberg correction applied for multiple testing
# for symmetric matrices, raw probabilites are reported below the diagonal and correlations adjusted for multiple comparisons above the diagonal
d5cordfbh <- corr.test(d5df, use ="pairwise", method="spearman", adjust="BH", alpha = .05)


###########################
## 16M-Delta-Correlation ##
###########################

## Subset all relevant baseline outcome measurements at t=6M
# The following primary endpoint variables are included in the analysis:
# Primary outcome variable: DM1ActivCV5
# Secondary outcome variables: 6MWTV5, PreBORGV5, AEScScoreV5, FDSSV5, MDHIV5, CISFatigueV5, BDIFsV5, INQOLQolScoreV5, stroopV5, MeanENMOV5, M5ENMOV5, L5ENMOV5
# Effect modifiers: MIRSV5, TMTV5, McGIllPainV5, ASBQV5, SSLDScoreV5, SSLNScoreV5, SSLIScoreV5, JFCSV5, ICQV5,
# IMQV5, SES28V5, CSIV5, AESIV5, CISactivityV5

# These variables are added seperately to adjust their names individually
dINQoLV5 <- df$INQOLQolScoreV5

V5df <- data.frame(
  df$DM1ActivCV5, 
  df$SMWTV5, df$PreBORGV5, df$AEScScoreV5, df$FDSSV5, df$MDHIV5, df$CISFatigueV5, df$BDIFsV5, dINQoLV5, dStroopInterferenceV5, df$MeanENMOV5, df$M5ENMOV5, df$L5ENMOV5,
  df$MIRSV5, dTMTV5, df$McGillPainV5, df$ASBQV5, df$SSLDScoreV5, df$SSLNScoreV5, df$SSLIScoreV5, df$JFCSV5, df$ICQV5,
  df$IMQV5, df$SES28V5, df$CSIV5, df$AESIV5, df$CISactivityV5)

# Check datatype of each variable
sapply(V5df, class) # -> MDHI, JFCS, ICQ, IMQ, SES28 loaded as 'factor'
V5df <- as.data.frame(sapply(V5df[1:length(colnames(V5df))], function(x){as.numeric(as.character(x))})) #changes all variable types to numeric

## Calculate delta-dataframe by subtracting bdf from edf
d16df <- V5df - bdf

## Clean up variable names
colnames(d16df) <- gsub(x = colnames(d16df), pattern = "df.", replacement = "d16")
colnames(d16df) <- gsub(x = colnames(d16df), pattern = "V5", replacement = "")
colnames(d16df) <- gsub(x = colnames(d16df), pattern = "d16DM1ActivC", replacement = "d16DM1-Activ-c")
colnames(d16df) <- gsub(x = colnames(d16df), pattern = "d16SMWT", replacement = "d166MWT")

## Restrict corrleation analysis to intervention group
d16df <- d16df[df$Treatment == "Intervention",] #n=128

## Check for missing values
sapply(d16df, function(x)sum(is.na(x)))

## Dropping variables with (almost) all values missing
d16df$d16SSLDScore <- NULL
d16df$d16SSLNScore <- NULL
d16df$d16SSLIScore <- NULL
d16df$d16JFCS <- NULL
d16df$d16ICQ <- NULL
d16df$d16IMQ <- NULL
d16df$d16SES28 <- NULL

## Shapiro-Wilks test for normality
pvalues <- apply(d16df, 2, function(x) shapiro.test(x)$p.value)
round(pvalues[pvalues > 0.05], 2) # 3/20 normally distributed

## Screening for outliers
# Outliers defined as < Q1-2*IQR | > Q3+2*IQR

cdf <- d16df #a copy is made in order to detect differnces after removal
for (x in 1:length(colnames(d16df))){
  Q1 <- unname(quantile(d16df[,x], 0.25, na.rm = T))
  Q3 <- unname(quantile(d16df[,x], 0.75, na.rm = T))
  IQR <- Q3-Q1
  l <- Q1 - 2*IQR
  u <- Q3 + 2*IQR
  d16df[,x][d16df[,x] > u | d16df[,x] < l] <- NA
}

# Number of values changed per column (=differences in NA values)
apply(cdf, 2, function(x) length(which(!is.na(x)))) - apply(d16df, 2, function(x) length(which(!is.na(x))))
sum(apply(cdf, 2, function(x) length(which(!is.na(x)))) - apply(d16df, 2, function(x) length(which(!is.na(x)))))
# -> 65 potential outliers identified, 12 removed

cdf$d166MWT[cdf$d166MWT < -100 | cdf$d166MWT > 210] <- NA #2
cdf$d16FDSS[cdf$d16FDSS > 30] <- NA #2
cdf$d16MDHI[cdf$d16MDHI < -50 | cdf$d16MDHI > 40] <- NA #2
cdf$d16BDIFs[cdf$d16BDIFs < -8] <- NA #2
cdf$d16ASBQ[cdf$d16ASBQ > 49] <- NA #2
cdf$d16CSI[cdf$d16CSI < - 5] <- NA #1
cdf$d16AESI[cdf$d16AESI < -15 | cdf$d16AESI > 20] <- NA #2

## restoring dataset and removing 12 outliers
d16df <- cdf

## Correlation analysis
# psych package using corr.test function
# non-parametric spearman correlations calculated given the non-normally distributed data
# Benjamini Hochberg correction applied for multiple testing
# for symmetric matrices, raw probabilites are reported below the diagonal and correlations adjusted for multiple comparisons above the diagonal
d16cordfbh <- corr.test(d16df, use ="pairwise", method="spearman", adjust="BH", alpha = .05)

###########################################
## Visualization of corr d5 and corr d16 ##
###########################################
plot.new()
par(mfrow=c(1,2))
corrplot(d5cordfbh$r, tl.col = "black", tl.srt = 45, 
         p.mat = d5cordfbh$p, sig.level = 0.05, insig = "blank")
mtext("A", side = 1, line = -45, cex=2, adj=0.05, font=2)

corrplot(d16cordfbh$r, tl.col = "black", tl.srt = 45, 
         p.mat = d16cordfbh$p, sig.level = 0.05, insig = "blank")
mtext("B", side = 1, line = -45, cex=2, adj=0.09, font=2)

# Exported using RStudio as 1900*1200 JPEG

##########################################
## Scatter plots of DM1-Activ-c vs 6MWT ##
##########################################

## Scatterplot of d5DM1ActivC and d56MWT
d5plot <- ggplot(d5df, aes_string(x=d5df$`d5DM1-Activ-c`, y=d5df$d56MWT)) +
  geom_point(color=("black"), size=2) +
  labs(x="d5DM1-Activ-c", y="d56MWT", tag="A") +
  scale_x_continuous(expand = c(0,0), limits = c(-25, 36), breaks = scales::pretty_breaks(n=10)) +
  scale_y_continuous(expand = c(0,0), limits = c(-90, 150), breaks = scales::pretty_breaks(n=10)) +
  geom_smooth(method="lm", color=("black"), se=F) +
  theme(
    panel.border = element_rect(colour="black", fill = NA, size = 0.4),
    axis.text = element_text(size = 15),
    axis.title = element_text(colour="black", size = 20, face="bold"),
    plot.tag.position = "topleft",
    plot.tag = element_text(size=30, face="bold"),
    panel.background = element_rect(fill="white"),
    aspect.ratio = 1
  )

## Scatterplot of d16DM1ActivC and d166MWT
d16plot <- ggplot(d16df, aes_string(x=d16df$`d16DM1-Activ-c`, y=d16df$d166MWT)) +
  geom_point(color=("black"), size=2) +
  labs(x="d16DM1-Activ-c", y="d166MWT", tag="B") +
  scale_x_continuous(expand = c(0,0), limits = c(-20, 36), breaks = scales::pretty_breaks(n=10)) +
  scale_y_continuous(expand = c(0,0), limits = c(-80, 210), breaks = scales::pretty_breaks(n=10)) +
  geom_smooth(method="lm", color=("black"), se=F) +
  theme(
    panel.border = element_rect(colour="black", fill = NA, size = 0.4),
    axis.text = element_text(size = 15),
    axis.title = element_text(colour="black", size = 20, face="bold"),
    plot.tag.position = "topleft",
    plot.tag = element_text(size=30, face="bold"),
    panel.background = element_rect(fill="white"),
    aspect.ratio = 1
  )
grid.arrange(d5plot, d16plot, ncol=2)

# Exported using RStudio as 1900*1200 JPEG









