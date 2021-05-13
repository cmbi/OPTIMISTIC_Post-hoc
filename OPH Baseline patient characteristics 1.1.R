# Script run with R version 4.0.5 (2021-03-31)
# Last changes applied on: 29/4/2020

###############
## Libraries ##
###############
library("readxl")

############################################
## Set working directory and open databse ##
############################################
setwd("Z:/documents")
df <- read_excel("Outcome measures overview (V5)_corrected.xlsx", sheet = 1)

## Location distribution
table(df$Centre)

# A: Paris: 71: 28%
# B: Munich: 66: 26%
# C: Nijmegen: 66: 26%
# D: Newcastle: 52: 20%

## Average age
mean(df$AgeBaseline)
# 45.6 years

## Onset age
mean(as.numeric(df$AgeAtOnset), na.rm = T)
# 25.4 years

## Sex distribution
table(df$Sex)
# 118 Female 46%
# 137 Male 54%

## CTG repeat length at baseline
mean(as.numeric(df$V2Mode), na.rm =T)
# 511

## Estimated CTG repeat length
mean(as.numeric(df$ePAL), na.rm =T)
# 257

