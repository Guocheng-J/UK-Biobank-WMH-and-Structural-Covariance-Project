###############################################################################
# Project description:
# This R script contains the following part of the code used in the manuscript
# " G Jiang et al. Common vascular white matter lesions contribute to specific 
#   grey matter network disruptions"
#
# Part I. Structural covariance network (SCN) simulations:
#    1. Creation of a 3x3 SCN adjacency matrix to inspect the step functions.
#        [Supplementary Figure 3 - Panel A]
#        You can plug in your own number and test it
#    2. Simulations: Assessment of effect size sensitivity for SCN connectivity 
#        comparisons.
#        [Supplementary Figure 5]
#    3. Simulation: Optimum folds of data for SCN analysis
#        [Supplementary Figure 4]
#
# Author: Guocheng Jiang
# Version: 2.1 (Last edits: Apr 24th 2026)
# R-version: 4.5.0
###############################################################################
# Loading required libraries
# To install missing R libraries: install.packages()
library(corrplot)          # To create visual for adjacency matrix
library(ggplot2)           # To create figures such as scatter plots
library(sn)                # To generate skew-normal distribution
library(survival)          # To create visual for SCN comparison
library(survminer)         # To create visual for SCN comparison

###############################################################################
# Project 1: 3x3 SCN simulation and step function visual inspection
###############################################################################
### 1.1 Function Definition ###

# 1.1a: Function to simulate symmetric adjacency matrix with optional skewness
#       Input: Matrix size, Mean and SD of covariance value, skewness
#       Skewness: 0 - Normal, Pos - right skew, Neg - left skew
#       e.g. Set a positively-skewed skewness values to simulate many weak 
#            connections with a few strong connections.

simulate_adjacency_matrix <- function(n, mean_val, sd_val, skewness = 0,
                                      min_val = -1, max_val = 1) {
    if (skewness == 0) {
    # Normal distribution
    values <- rnorm(n * n, mean = mean_val, sd = sd_val)
    
  } else {
    # Skewed distribution using rsn function
    values <- rsn(n * n, xi = mean_val, omega = sd_val, alpha = skewness)
    
  }
  mat <- matrix(values, nrow = n, ncol = n)
  
  # Make the matrix symmetric
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  mat <- pmax(pmin(mat, max_val), min_val)
  # Remove self-connections
  diag(mat) <- NA                                                                
  return(mat)
}

# 1.1b: Function to extract values from adjacency matrix and generate data for 
#       step function plotting. Here we used a survival plot library to help 
#       visualization
#       Input: Adjacency Matrix, Threshold of matrix index
#       Output: A "survival" data frame with threshold (as "time") and removal of 
#               the edge < threshold (as "event")

extract_thresholded_SCN_data <- function(matrix, thresholds) {
  # Extract only lower triangle values since it is symmetric.
  values <- na.omit(matrix[lower.tri(matrix)])
  thresholded_SCN_data <- data.frame()
  
  # Determine the step differences in the threshold
  dt   <- diff(thresholds)[1] 
  
  # Create a data frame to store the thresholded SCN data:
  for (thresh in thresholds) {
    # Define an edge survives if its value > threshold
    survived <- sum(values > thresh)
    # Define an edge die if its value < threshold, and count new "dead" edges
    deaths <- sum(values <= thresh & values > (thresh - dt)) 
    thresholded_SCN_data <- rbind(thresholded_SCN_data, data.frame(Threshold = thresh, Survived = survived, Deaths = deaths))
  }
  
  return(thresholded_SCN_data)
}

# 1.1c: Function to format data for Kaplan-Meier plots with censor handles:
#       Input: The thresholded SCN dataframe, and a label of groups

format_km_data <- function(thresholded_SCN_data, group_label) {
  df <- data.frame(
    # Expand based on number of deaths
    time = rep(thresholded_SCN_data$Threshold, thresholded_SCN_data$Deaths), 
    # Count each death is an event
    status = rep(1, sum(thresholded_SCN_data$Deaths)),
    # Creating group labels
    group = rep(group_label, sum(thresholded_SCN_data$Deaths))
  )
  
  # Handle censoring: If an edge survives beyond the last threshold, add it as censored
  max_threshold <- max(thresholded_SCN_data$Threshold)
  # Last number of survived edges
  final_survived <- thresholded_SCN_data$Survived[nrow(thresholded_SCN_data)]
  if (final_survived > 0) {
    censoring_df <- data.frame(
      # Repeat for each surviving edge
      time = rep(max_threshold, final_survived),
      # Update the status to "Censored". All edges survived the max threshold
      # (in manuscript: 0.5) is treated as censored data.
      status = rep(0, final_survived),
      group = rep(group_label, final_survived)
    )
    df <- rbind(df, censoring_df)
  }
  
  return(df)
}


#------------------------------------------------------------------------------
### 1.2 Simulation of one 3x3 adjacency matrix ###
# You can put down any number between -1 and 1 to create your own SCN.

# SCN adjacency matrix for Group A (e.g. A disrupted SCN)
test_mat1 <- matrix(c(
  1.0,  0.5, 0.01,   # correlations for Var1
  0.25, 1.0,  0.3,   # correlations for Var2
  0.15, 0.01, 1.0    # correlations for Var3
), nrow = 3, byrow = TRUE)

# Rename your rows/columns
colnames(test_mat1) <- rownames(test_mat1) <- c("Var1", "Var2", "Var3")

# SCN adjacency matrix for Group B (e.g. A well-preserved SCN)
test_mat2 <- matrix(c(
  1.0,  0.5, -0.2,   # correlations for Var1
  0.5,  1.0,  0.3,   # correlations for Var2
  0.45, 0.4,  1.0    # correlations for Var3
), nrow = 3, byrow = TRUE)

# Rename your rows/columns
colnames(test_mat2) <- rownames(test_mat2) <- c("Var1", "Var2", "Var3")

# Remove the diagnal (self-correlation) elements
diag(test_mat2) <- 0
diag(test_mat1) <- 0

# Visualization of the adjacency matrix:
custom_colors <- colorRampPalette(c("navy","royalblue3", "white", "firebrick1", "black"))(200)  # Deep color version
corrplot(test_mat1, method = "color", type = "lower", col = custom_colors, addCoef.col = "black", tl.col = "black", tl.cex = 1)
corrplot(test_mat2, method = "color", type = "lower", col = custom_colors, addCoef.col = "black", tl.col = "black", tl.cex = 1)

# Extract the K-M matrix
thresh_range <- seq(0, 0.5, by = 0.01)  # Threshold value
Mat1_thres <- extract_thresholded_SCN_data(test_mat1, thresh_range) # Mat1 have many low correlated volumes
Mat2_thres <- extract_thresholded_SCN_data(test_mat2, thresh_range) # Mat2 have many high correlated volumes
Mat1_km <- format_km_data(Mat1_thres, "Matrix1_sim")
Mat2_km <- format_km_data(Mat2_thres, "Matrix2_sim")
simulation_km_data <- rbind(Mat1_km, Mat2_km)

# Plot the step function using Kaplan Meier method from ggsurvplot function
km_fit <- survfit(Surv(time, status) ~ group, data = simulation_km_data)
ggsurvplot(km_fit, data = simulation_km_data, pval = FALSE,
           legend.title = "Group",
           xlab = "Pearson Correlation Threshold",
           ylab = "Perecent Connections Remaining",
           conf.int = FALSE,
           risk.table = FALSE,
           size = 2,
           ggtheme = theme_bw())




###############################################################################
# Project 2: Extend the idea to NxN adjacency matrix (SCN with N nodes)
#            and examine the minimal effect size needed for stats model to 
#            reveal significant differences between the two SCNs. 
#            (e.g. Log rank test, and Cox Proportional Hazards Model)
###############################################################################
# Simulation parameters
# You can adjust the parameters to run your own simulation.
# Note: This step will be computationally intensive and will take hours 
#       depending on how many simulations you will perform.
thresh_range <- seq(0, 0.5, by = 0.01)         # Range of threshold
n_values <- seq(8, 50, by = 1)                 # Size of SCN (num. of nodes)
num_simulations <- 100                         # Number of random iterations
power_threshold <- 0.8                         # Type 2 error rate
ctrl_effect <- 0.2                             # SCN average weights at Control Group
sd_value <- 0.1                                # SD of the matrix elements
skew_value <- 10                               # Skewness

# Create a data frame to store results
results_power_analysis <- data.frame(N = integer(), Min_Cohen_D = numeric())

# Before running the whole simulation, you can benchmark yourself by creating some 
# random matrix to examine your parameter input

n_test <- 100
test_matrix <- simulate_adjacency_matrix(n_test, mean_val = ctrl_effect, sd_val = sd_value, skewness = skew_value)
custom_colors <- colorRampPalette(c("navy","royalblue3", "white", "firebrick1", "black"))(200)  # Deep color version
corrplot(test_matrix, method = "color", type = "lower", col = custom_colors, tl.col = "black", tl.cex = 1)
hist(test_matrix)   # The histogram shows the distribution of edge weights

#------------------------------------------------------------------------------
# 2.2 Simulations
# Loop over different ROI sizes
for (n in n_values) {
  cat("Processing N =", n, "\n")
  cohen_d_values <- c()
  
  # Running n numbers of simulations
  for (replicate in 1:num_simulations) {
    cat("Replication =", replicate , "\n")
    for (effect_size in seq(0, 0.3, by = 0.01)) {
      cat("    >> Simulating effect size =", effect_size , "\n")
      significant_count <- 0
      
      for (i in 1:num_simulations) {
        # Generate matrices
        control_matrix <- simulate_adjacency_matrix(n, mean_val = ctrl_effect, sd_val = sd_value, skewness = skew_value)
        disease_matrix <- simulate_adjacency_matrix(n, mean_val = ctrl_effect - effect_size, sd_val = sd_value, skewness = skew_value)
        
        # Extract survival data
        control_survival <- extract_thresholded_SCN_data(control_matrix, thresh_range)
        disease_survival <- extract_thresholded_SCN_data(disease_matrix, thresh_range)
        
        # Format for survival analysis
        control_km <- format_km_data(control_survival, "Control")
        disease_km <- format_km_data(disease_survival, "Disease")
        
        # Combine datasets
        survival_km_data <- rbind(control_km, disease_km)
        
        # Conduct Log-Rank Test, you can swap your own stats tests here
        # e.g. Cox model, don't forget redefine p-value extractions
        logrank_test <- survdiff(Surv(time, status) ~ group, data = survival_km_data)
        p_value <- 1 - pchisq(logrank_test$chisq, df = 1)
        
        # Count significant results
        if (p_value < 0.05) {
          significant_count <- significant_count + 1
        }
      }
      
      # Compute power
      power <- significant_count / num_simulations
      
      # If power > 0.8, compute Cohen's D
      if (power >= power_threshold) {
        mean_control <- mean(control_matrix[lower.tri(control_matrix)], na.rm = TRUE)
        mean_disease <- mean(disease_matrix[lower.tri(disease_matrix)], na.rm = TRUE)
        pooled_sd <- sqrt((var(control_matrix[lower.tri(control_matrix)], na.rm = TRUE) + 
                             var(disease_matrix[lower.tri(disease_matrix)], na.rm = TRUE)) / 2)
        cohen_d <- (mean_control - mean_disease) / pooled_sd
        cohen_d_values <- c(cohen_d_values, cohen_d)
        break
      }
    }
  }
  
  # Compute confidence intervals
  if (length(cohen_d_values) > 0) {
    min_cohen_d <- mean(cohen_d_values)
    ci_lower <- quantile(cohen_d_values, 0.025)
    ci_upper <- quantile(cohen_d_values, 0.975)
    
    results_power_analysis <- rbind(results_power_analysis, 
                                    data.frame(N = n, Min_Cohen_D = min_cohen_d, CI_Lower = ci_lower, CI_Upper = ci_upper))
  }
}

#------------------------------------------------------------------------------
# 2.3 Visualization of the results
# Plot results with confidence intervals
ggplot(results_power_analysis, aes(x = N, y = Min_Cohen_D)) +
  geom_point(size = 2, color = "red") +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), fill = "lightblue", alpha = 0.3) +
  labs(title = "Minimum Cohen's D for >80% Power to show significant differences in SCNs",
       x = "Number of ROIs (Matrix Size)",
       y = "Minimum Cohen's D for 80% Power") +
  theme_bw()


###############################################################################
# Project 3: Effect of data stratification on the estimation of age-related 
#            differences in SCN connectivity. Effects on global SCN 
#            connectivity between the sixties and seventies groups were 
#            estimated using Cohen’s D. 
#
# Note: You will need to plug in your own data. In the manuscript we examined
#       the age effect on SCN. This simulator allows you to create data folds,
#       and examine the optimum folds you can stratified your data into to 
#       quantify the within-group variabilities in the SCN connectivity
#       
###############################################################################
#------------------------------------------------------------------------------
# 3.1 Function definition:
#
# 3.1a To create a dataframe of normalized volumetric data. The data frame will
#      be used to generate adjacency matrix for SCN analysis. 
#
#      Dependent variable names need to match the names of ROIs on your datasheet.
#     To create your customized SCN graphs, you can add/drop nodes in this function.
#      
#
# In this script, the following independent and confounding variables were used:
#   Age_baseline_imaging: Chronological age (years) at baseline scan.
#   Sex_value: Biological sex categorical variable. 0=Female, 1=Male.
#   eTIV_baseline_L: Estimated intracranial volume (Unit: L)
#   vCSF_vol_ml: Ventricular cerebrospinal fluid volume (Unit: ml)
#
#
GenerateNode <- function(Graph_input) {
  Nodes <- data.frame(
    Amyg = lm(Amyg_v ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    AG = lm(Angular_GMV_bilateral ~ Age_baseline_imaging + 
              Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
            data=Graph_input)$residuals,
    Caud = lm(Caud_v_1 ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    COC = lm(Central_opercular_cortex_GMV_bilateral ~ Age_baseline_imaging + 
               Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
             data=Graph_input)$residuals,
    aCG = lm(ACC_GMV_bilateral ~ Age_baseline_imaging + 
               Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
             data=Graph_input)$residuals,
    pCG = lm(PCC_GMV_bilateral ~ Age_baseline_imaging + 
               Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
             data=Graph_input)$residuals,
    FMC =  lm(Frontal_medial_cortex_GMV_bilateral ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    FOpe =  lm(Frontal_operculum_GMV_bilateral ~ Age_baseline_imaging + 
                 Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
               data=Graph_input)$residuals,
    OFC =  lm(Frontal_Orbital_GMV_bilateral ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    FP =  lm(Frontal_Pole_GMV_bilateral ~ Age_baseline_imaging + 
               Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
             data=Graph_input)$residuals,
    Hipp =  lm(Hipp_v ~ Age_baseline_imaging + 
                 Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
               data=Graph_input)$residuals,
    IFPO =  lm(Inferior_Frontal_opercularis_GMV_bilateral ~ Age_baseline_imaging + 
                 Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
               data=Graph_input)$residuals,
    IFPT = lm(Inferior_Frontal_triangularis_GMV_bilateral ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    Insu = lm(Insula_GMV_bilateral ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    MFG = lm(Middle_Frontal_GMV_bilateral ~ Age_baseline_imaging + 
               Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
             data=Graph_input)$residuals,
    Pall = lm(Pall_v ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    PaCG = lm(Paracingulate_GMV_bilateral ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    aPHG = lm(ParaHippo_ant_GMV_bilateral ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    pPHG = lm(ParaHippo_post_GMV_bilateral ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    POC = lm(Parietal_operculum_GMV_bilateral ~ Age_baseline_imaging + 
               Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
             data=Graph_input)$residuals,
    PP = lm(Planum_Polare_bilateral ~ Age_baseline_imaging + 
              Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
            data=Graph_input)$residuals,
    Prec = lm(Precuneous_GMV_bilateral ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    Puta = lm(Puta_v ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    Scal = lm(Subcallosal_GMV_bilateral ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    SFG = lm(Superior_Frontal_GMV_bilateral ~ Age_baseline_imaging + 
               Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
             data=Graph_input)$residuals,
    aSMG = lm(Supramarginal_ant_GMV_bilateral ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    pSMG = lm(Supramarginal_post_GMV_bilateral ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    Thal = lm(Thal_v_1 ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    NAcc = lm(NAcc_v ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    aITG = lm(Inferior_Temporal_Ant_GMV_bilateral ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    pITG = lm(Inferior_Temporal_Post_GMV_bilateral ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    toITG = lm(Inferior_Temporal_Tempocc_GMV_bilateral ~ Age_baseline_imaging + 
                 Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
               data=Graph_input)$residuals,
    aSTG = lm(Superior_Temporal_ant_GMV_bilateral ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    pSTG = lm(Superior_Temporal_post_GMV_bilateral ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    aMTG = lm(Middle_Temporal_ant_GMV_bilateral ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    pMTG = lm(Middle_Temporal_post_GMV_bilateral ~ Age_baseline_imaging + 
                Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    toMTG = lm(Middle_Temporal_tempocc_GMV_bilateral ~ Age_baseline_imaging + 
                 Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
               data=Graph_input)$residuals
    
  )
  return(Nodes)
}

# 3.1b The same code but removed age covariate in the study to demonstrate age effect.
GenerateNode_AgeEffect <- function(Graph_input) {
  Nodes <- data.frame(
    Amyg = lm(Amyg_v ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    AG = lm(Angular_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
            data=Graph_input)$residuals,
    Caud = lm(Caud_v_1 ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    COC = lm(Central_opercular_cortex_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
             data=Graph_input)$residuals,
    aCG = lm(ACC_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
             data=Graph_input)$residuals,
    pCG = lm(PCC_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
             data=Graph_input)$residuals,
    FMC =  lm(Frontal_medial_cortex_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    FOpe =  lm(Frontal_operculum_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
               data=Graph_input)$residuals,
    OFC =  lm(Frontal_Orbital_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    FP =  lm(Frontal_Pole_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
             data=Graph_input)$residuals,
    Hipp =  lm(Hipp_v ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
               data=Graph_input)$residuals,
    IFPO =  lm(Inferior_Frontal_opercularis_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
               data=Graph_input)$residuals,
    IFPT = lm(Inferior_Frontal_triangularis_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    Insu = lm(Insula_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    MFG = lm(Middle_Frontal_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
             data=Graph_input)$residuals,
    Pall = lm(Pall_v ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    PaCG = lm(Paracingulate_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    aPHG = lm(ParaHippo_ant_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    pPHG = lm(ParaHippo_post_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    POC = lm(Parietal_operculum_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
             data=Graph_input)$residuals,
    PP = lm(Planum_Polare_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
            data=Graph_input)$residuals,
    Prec = lm(Precuneous_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    Puta = lm(Puta_v ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    Scal = lm(Subcallosal_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    SFG = lm(Superior_Frontal_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
             data=Graph_input)$residuals,
    aSMG = lm(Supramarginal_ant_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    pSMG = lm(Supramarginal_post_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    Thal = lm(Thal_v_1 ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    NAcc = lm(NAcc_v ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    aITG = lm(Inferior_Temporal_Ant_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    pITG = lm(Inferior_Temporal_Post_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    toITG = lm(Inferior_Temporal_Tempocc_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
               data=Graph_input)$residuals,
    aSTG = lm(Superior_Temporal_ant_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    pSTG = lm(Superior_Temporal_post_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    aMTG = lm(Middle_Temporal_ant_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    pMTG = lm(Middle_Temporal_post_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
              data=Graph_input)$residuals,
    toMTG = lm(Middle_Temporal_tempocc_GMV_bilateral ~ Sex_value + eTIV_baseline_L + vCSF_vol_ml, 
               data=Graph_input)$residuals
    
  )
  return(Nodes)
}

# 3.1c Function to compute adjacency matrix
Compute_AdjacencyMatrix <- function(data) {
  cor_matrix <- cor(data)                             # Calculate the correlation matrix
  adjacency_matrix <- cor_matrix                      # Create an adjacency matrix
  diag(adjacency_matrix) <- 0                         # Set diagonal to 0 (no self-loops)
  return(adjacency_matrix)
}

# 3.1d Effect size estimation from adjacency matrix
cohen_d_from_adj <- function(adj_young, adj_old) {
  # Use lower triangle edges as "observations"
  e_y <- adj_young[lower.tri(adj_young)]
  e_o <- adj_old[lower.tri(adj_old)]
  
  # Remove NA if any
  e_y <- e_y[is.finite(e_y)]
  e_o <- e_o[is.finite(e_o)]
  
  if (length(e_y) < 5 || length(e_o) < 5) return(NA_real_)
  
  m_y <- mean(e_y)
  m_o <- mean(e_o)
  v_y <- var(e_y)
  v_o <- var(e_o)
  
  pooled_sd <- sqrt((v_y + v_o) / 2)
  if (!is.finite(pooled_sd) || pooled_sd == 0) return(NA_real_)
  
  (m_o - m_y) / pooled_sd  # positive => Older higher mean edge than Younger
}


#------------------------------------------------------------------------------
# 3.2 Data input
# Here is an example of input sixties vs. seventies in our data sheet.
# In our age effect study, we remove the confounder "Age_baseline_imaging" in
# the GenerateNode() function.

Node_Sixties <- GenerateNode_AgeEffect(UKB_T1wImaging_Vol[UKB_T1wImaging_Vol$AgeGroup == "Sixties", ])
Node_Seventies <- GenerateNode_AgeEffect(UKB_T1wImaging_Vol[UKB_T1wImaging_Vol$AgeGroup == "Seventies", ])

# Define the number of data folds you want to simulate: e.g. 2, 4, 8, 16 folds.
K_values <- c(2, 4, 8, 16, 32)


#------------------------------------------------------------------------------
# 3.3 SCN construction and effect size calculation at full sample
# Create adjacency matrices
AdjMat_Sixties <- Compute_AdjacencyMatrix(Node_Sixties)
AdjMat_Seventies <- Compute_AdjacencyMatrix(Node_Seventies)

# Calculate the effect size when supporting only 1 fold of the data:
cohen_d_from_adj(AdjMat_Sixties, AdjMat_Seventies)


#------------------------------------------------------------------------------
# 3.4 Estimation of effect size distribution at different folds of data
# Split out Younger and Older once
df_sixties <- UKB_T1wImaging_Vol %>% filter(AgeGroup == "Sixties")
df_seventies <- UKB_T1wImaging_Vol %>% filter(AgeGroup == "Seventies")


# Initiate a result dataframe as a list.
results <- list()

# Loop through different k values:
for (K in K_values) {
  cat("Processing K =", K, "\n")
  
  # 1) Create random folds *within each group* (balanced)
  fold_y <- sample(rep(1:K, length.out = nrow(df_sixties)))
  fold_o <- sample(rep(1:K, length.out = nrow(df_seventies)))
  
  # 2) For each fold, compute SCN matrices and Cohen's D
  d_vec <- sapply(1:K, function(k) {
    df_yk <- df_sixties[fold_y == k, , drop = FALSE]
    df_ok <- df_seventies[fold_o == k, , drop = FALSE]
    
    # Optional sanity check: avoid tiny folds (can remove if you want)
    if (nrow(df_yk) < 5 || nrow(df_ok) < 5) return(NA_real_)
    
    Node_Y <- GenerateNode_AgeEffect(df_yk)
    Node_O <- GenerateNode_AgeEffect(df_ok)
    
    Adj_Y <- Compute_AdjacencyMatrix(Node_Y)
    Adj_O <- Compute_AdjacencyMatrix(Node_O)
    
    cohen_d_from_adj(Adj_Y, Adj_O)
  })
  
  results[[as.character(K)]] <- data.frame(
    K = K,
    fold = 1:K,
    cohen_d = d_vec
  )
}

# Transfer the simulation results into the output dataframe
results_df <- bind_rows(results)

# Plot distribution of fold-wise Cohen's D by K
ggplot(results_df, aes(x = factor(K), y = cohen_d)) +
  #geom_violin(trim = TRUE) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  labs(
    title = "Instability of fold-wise Cohen's D with increasing folds",
    x = "Number of folds (K)",
    y = "Cohen's D (edge-distribution; Sixties vs Seventies group)"
  ) +
  theme_bw()




