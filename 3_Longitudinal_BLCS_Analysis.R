###############################################################################
# Project description:
# This R script contains the following part of the code used in the manuscript
# " G Jiang et al. Common vascular white matter lesions contribute to specific 
#   grey matter network disruptions"
#
# Part III. Longitudinal Analysis using Bivariate Latent Change Score Models
#
#
# Author: Guocheng Jiang
# Version: 2.1 (Last edits: Apr 24th 2026)
# R-version: 4.5.0
###############################################################################
# Loading required libraries
# To install missing R libraries: install.packages()
library(lavaan)
library(dplyr)

#------------------------------------------------------------------------------
# 1. Create the longitudinal data frame:
# In our manuscript, we named the tabular data as "older". Here you will need to
# feed in your own data frame and adjust the column names. In this example, we 
# demonstrated method to conduct BLCS analysis on the bidirectional relationship
# between white matter hyperintensity (WMH) accumulation and grey matter volume
# loss in insular cortex.

summary(older$Insula_GMV_bilateral_2)
older$Insula_GMV_bilateral_1 <- older$Insula_GMV_bilateral /1000
older$Insula_GMV_bilateral_2 <- older$Insula_GMV_bilateral_2 /1000

Longitudinal_BLCS_frame <- older %>%
  filter(
    # Primary variables of interests
    !is.na(Insula_GMV_bilateral_2),
    !is.na(Insula_GMV_bilateral_1),
    !is.na(WMH_vol_baseline_log_ml),
    !is.na(WMH_vol_followup_log_ml),
   
    # Confounding variables
    !is.na(Age_baseline_imaging),
    !is.na(Sex_value),
    !is.na(eTIV_baseline_L)
  )


# Anatomical volume normalization in both baseline and repeat imaging estimates
# Regression with confounding variables including age, sex, eTIV, and vCSF volume. 
Longitudinal_BLCS_frame$Insula_GMV_bilateral_1_normalized = lm(Longitudinal_BLCS_frame$Insula_GMV_bilateral_1 ~ 
                                                                 Longitudinal_BLCS_frame$Age_baseline_imaging +
                                                                 Longitudinal_BLCS_frame$Sex_value +
                                                                 Longitudinal_BLCS_frame$eTIV_baseline_L)$residuals

Longitudinal_BLCS_frame$Insula_GMV_bilateral_2_normalized = lm(Longitudinal_BLCS_frame$Insula_GMV_bilateral_2 ~ 
                                                                 Longitudinal_BLCS_frame$Age_baseline_imaging +
                                                                 Longitudinal_BLCS_frame$Sex_value +
                                                                 Longitudinal_BLCS_frame$eTIV_baseline_L)$residuals

Longitudinal_BLCS_frame$WMH_vol_baseline_log_1_normalized = lm(Longitudinal_BLCS_frame$WMH_vol_baseline_log_ml ~ 
                                                                 Longitudinal_BLCS_frame$Age_baseline_imaging +
                                                                 Longitudinal_BLCS_frame$Sex_value +
                                                                 Longitudinal_BLCS_frame$eTIV_baseline_L)$residuals

Longitudinal_BLCS_frame$WMH_vol_baseline_log_2_normalized = lm(Longitudinal_BLCS_frame$WMH_vol_followup_log_ml ~ 
                                                                 Longitudinal_BLCS_frame$Age_baseline_imaging +
                                                                 Longitudinal_BLCS_frame$Sex_value +
                                                                 Longitudinal_BLCS_frame$eTIV_baseline_L)$residuals

# Prepare BLCS entry
Ourdata <- dplyr::select(Longitudinal_BLCS_frame, Insula_GMV_bilateral_2_normalized, WMH_vol_baseline_log_2_normalized,
                                                  Insula_GMV_bilateral_1_normalized, WMH_vol_baseline_log_1_normalized)
colnames(Ourdata) <- c("InsulaVol2", "WMHVol2", "InsulaVol1", "WMHVol1")

#------------------------------------------------------------------------------
# 1. Create the longitudinal data frame:
# Entry of BLCS model:

BLCS<-'

InsulaVol2 ~ 1*InsulaVol1       # This parameter regresses InsulaVol2 perfectly on InsulaVol1
dInsulaVol1 =~ 1*InsulaVol2     # This defines the latent change score factor as measured perfectly by scores on InsulaVol2
dInsulaVol1 ~ 1                 # This estimates the intercept of the change score 
InsulaVol1 ~  1                 # This estimates the intercept of InsulaVol1 
InsulaVol2 ~ 0*1                # This constrains the intercept of InsulaVol2 to 0

WMHVol2 ~ 1*WMHVol1             # This parameter regresses WMHVol2 perfectly on WMHVol1
dWMH1 =~ 1*WMHVol2              # This defines the latent change score factor as measured perfectly by scores on WMHVol2
WMHVol2 ~ 0*1                   # This line constrains the intercept of WMHVol2 to 0
WMHVol2 ~~ 0*WMHVol2            # This fixes the variance of the WMHVol1 to 0  

dInsulaVol1 ~~  dInsulaVol1     # This estimates the variance of the change scores
InsulaVol1 ~~   InsulaVol1      # This estimates the variance of the InsulaVol1 
InsulaVol2 ~~ 0*InsulaVol2      # This fixes the variance of the InsulaVol2 to 0  

dWMH1 ~ 1                       # This estimates the intercept of the change score 
WMHVol1 ~ 1                     # This estimates the intercept of WMHVol1 
dWMH1 ~~ dWMH1                  # This estimates the variance of the change scores 
WMHVol1 ~~ WMHVol1              # This estimates the variance of WMHVol1 

dWMH1~InsulaVol1 + WMHVol1      # This estimates the InsulaVol to WMHVol coupling parameter and the InsulaVol to InsulaVol self-feedback
dInsulaVol1~WMHVol1 + InsulaVol1# This estimates the WMHVol to InsulaVol coupling parameter and the WMHVol to WMHVol self-feedback

InsulaVol1 ~~  WMHVol1          # This estimates the InsulaVol1 WMHVol1 covariance
dInsulaVol1~~dWMH1              # This estimates the dInsulaVol and dWMHVol covariance
'

# Model fitting: BLCS with propensity matched sample
fitBLCS <- lavaan(BLCS, data=Ourdata, estimator='mlr',fixed.x=FALSE,missing='fiml')
summary(fitBLCS, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE)

# You may import the results to open-access softwares such as Ωnyx for visualization.
# We acknowledge the software developments and related tutorials from the following
# works from Rogier A. Kievit et al. (2018)

# Citation: 
# Kievit, Rogier A, et al. “Developmental Cognitive Neuroscience Using Latent 
# Change Score Models: A Tutorial and Applications.” Developmental Cognitive 
# Neuroscience [Netherlands], vol. 33, Oct. 2018, pp. 99–117,
# https://doi.org/10.1016/j.dcn.2017.11.007.

