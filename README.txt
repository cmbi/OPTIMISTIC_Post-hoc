R-scripts used to generate the results as presented in the OPTIMISTIC post-hoc paper: "Clinical outcome evaluations and CBT response prediction in Myotonic Dystrophy" as well as non-published supplementary results

Scripts:
*OPH Baseline patient characteristics 1.1: Used for exploratory dataset analyses
*OPH Correlation analyses 1.2: Cross-sectional as well as intervention driven longitudinal Spearman's Rank correlations of outcome measurements used in the OPTIMISTIC trial.
*OPH Model validation - Bootstrap enhanced Elastic-Net 1.3: Reproduction of published DM1-Activ-c prediction results and comparison of models derived by Bootstrap enhanced Elastic-Net in different downsizing samples.

*OPH d10M_DM1Activc_prediction 1.1: Idenfitication of variables associated with changes in the DM1-Activ-c questionnare over the course of 10 months within the OPTIMISTIC intervention cohort using Bootstrap enhanced Elastic-Net regression analysis
*OPH d10M_6MWT_prediction 1.1: Idenfitication of variables associated with changes in 6MWT scores over the course of 10 months within the OPTIMISTIC intervention cohort using Bootstrap enhanced Elastic-Net regression analysis

*OPH d10M_DM1Activc_prediction - sensitivity control 1.1.R: Identification of variables associated with 10-month natural history (assessed by DM1-Activ-c) in the OPTIMISTIC control group 
*OPH d10M_6MWT_prediction - sensitivity control.R: Identification of variables associated with 10-month natural history (assessed by 6MWT) in the OPTIMISTIC control group

Results:
* OPH_DM1Activc_VIPs.csv: Variable Inclusion Probabilities of delta-DM1-Activ-c predictors among 10 imputed datasets [intervention group]
* OPH_DM1Activc_senscontrol_VIPs.csv: Variable Inclusion Probabilities of delta-DM1-Activ-c among 10 imputed datasets [control group]
* OPH_6MWT_VIPs.csv: Variable Inclusion Probabilities of delta-6MWT predictors among 10 imputed datasets [intervention group]
* OPH_6MWT_senscontrol_VIPs.csv: Variable Inclusion Probabilities of delta-6MWT predictors among 10 imputed datasets [control group]
* OPH_delta_6MWT_fit: Regression fit on unimputed intervention cases using the selected predictors with avg VIPs > 60%
