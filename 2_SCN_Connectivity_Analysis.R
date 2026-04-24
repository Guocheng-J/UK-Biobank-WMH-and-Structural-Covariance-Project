###############################################################################
# Project description:
# This R script contains the following part of the code used in the manuscript
# " G Jiang et al. Common vascular white matter lesions contribute to specific 
#   grey matter network disruptions"
#
# Part II. Structural covariance network (SCN) analysis: Comparisons of overall
#          and backbone network connectivity between two groups.
#
#       1. Create SCN adjacency matrix given a data frame of volumetric estimates
#          and compare global SCN connectivity between two groups.
# 
#       2. Estimate backbone network connections in the SCN graph and visualize
#          the backbone network connectivity difference using ribbon and 
#          circular diagrams.
#
# Author: Guocheng Jiang
# Version: 2.1 (Last edits: Apr 24th 2026)
# R-version: 4.5.0
###############################################################################
# Loading required libraries
# To install missing R libraries: install.packages()

library(corrplot)          # To create visual for adjacency matrix
library(circlize)          # To visualize backbone network organization
library(dplyr)             # To create data strata and perform tabular operation
library(ggplot2)           # To create figures such as scatter plots
library(igraph)            # To conduct graph theory calculations
library(tidyr)             # To perform data cleaning and reshaping
library(scales)            # To linearly maps variables into a new range
library(sn)                # To generate skew-normal distribution
library(stringr)           # To detect keywords in a string
library(survival)          # To create visual for SCN comparison
library(survminer)         # To create visual for SCN comparison

###############################################################################
# Section 1: Create SCN adjacency matrices from volumetric tabular data
#            Participants were ranked within each single-year age stratum 
#            (e.g. 65, 66, 67… years) based on the transformed and normalized 
#            WMH volume and were then stratified into WMH-low (<30 percentile) 
#            or WMH-high (>70 percentile) groups. 
#
#            The same code base was used to compare age effect, sex effect, and 
#            vCSF effect by replacing the input data and stratification 
#            criteria
###############################################################################
#------------------------------------------------------------------------------
# Function definitions:

# 1.   To create a dataframe of normalized volumetric data. The data frame will
#      be used to generate adjacency matrix for SCN analysis. 
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
    #Cune = lm(Cuneal_cortex_GMV_bilateral ~ Age_baseline_imaging + 
    #           Sex_value + eTIV_baseline_L, 
    #         data=Graph_input)$residuals,
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

# 2.   To compute adjacency matrix
Compute_AdjacencyMatrix <- function(data) {
  cor_matrix <- cor(data)                             # Calculate the correlation matrix
  adjacency_matrix <- cor_matrix                      # Create an adjacency matrix
  diag(adjacency_matrix) <- 0                         # Set diagonal to 0 (no self-loops)
  return(adjacency_matrix)
}

# 3.   Function to extract values from adjacency matrix and generate data for 
#      step function plotting. Here we used a survival plot library to help 
#      visualization.
#      Input: Adjacency Matrix, Threshold of matrix index
#      Output: A "survival" data frame with threshold (as "time") and removal of 
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



#------------------------------------------------------------------------------
# 1.1 Group stratification and abbreviation definition
# 

# Loading the data
older <- data %>% filter(Age_baseline_imaging >= 65)

# Stratification into WMH High and WMH Low Group

older <- older %>%
  # Participants were ranked within each single year age stratum:
  # i.e. 65 years old compare to 65 years old.
  group_by(Age_baseline_imaging) %>%
  mutate(
    # Calculate the 30th and 70th percentile of WMH for each age group
    p30 = quantile(WMH_log_vol_regressed, 0.3, na.rm = TRUE),
    p70 = quantile(WMH_log_vol_regressed, 0.7, na.rm = TRUE),
    
    # Assign WMH status based on percentiles
    `WMH_status` = case_when(
      WMH_log_vol_regressed > p70 ~ "High",
      WMH_log_vol_regressed < p30 ~ "Low",
      WMH_log_vol_regressed >= p30 & WMH_log_vol_regressed <= p70 ~ "Medium",
      TRUE ~ ""
    )
  ) %>%
  ungroup()  # Ungroup after applying the operations

#------------------------------------------------------------------------------
# 1.2 Create adjacency matrix
Node_WMH_High <- GenerateNode(older[older$WMH_status == "High", ])
Node_WMH_Low <- GenerateNode(older[older$WMH_status == "Low", ])
AdjMat_WMH_High <- Compute_AdjacencyMatrix(Node_WMH_High)
AdjMat_WMH_Low <- Compute_AdjacencyMatrix(Node_WMH_Low)

# Visually inspect the plot
custom_colors <- colorRampPalette(c("navy","royalblue3", "white", "firebrick1", "black"))(200) 
corrplot(AdjMat_WMH_High, method = "color", type = "lower", col = custom_colors, 
         tl.col = "black", tl.cex = 1)
corrplot(AdjMat_WMH_Low, method = "color", type = "lower", col = custom_colors, 
         tl.col = "black", tl.cex = 1)


#------------------------------------------------------------------------------
# 1.3 Calculate percent of SCN edges remained after each threshold
# Initialize data storage
all_survival_data <- data.frame()
all_neg_survival_data <- data.frame()
logrank_results <- data.frame(Fold = integer(), Chi_Square = numeric(), P_Value = numeric())
neg_logrank_results <- data.frame(Fold = integer(), Chi_Square = numeric(), P_Value = numeric())

# Define threshold range
thresh_range <- seq(0, 0.5, by = 0.01)

# Create stratification
WMH_HighGroup <- older %>% filter(WMH_status == "High")
WMH_LowGroup <- older %>% filter(WMH_status == "Low")

# Define number of folds
n_folds <- 8

# Assign random folds:
WMH_HighGroup <- WMH_HighGroup %>%
  mutate(Fold = sample(rep(1:n_folds, length.out = n()), n())) %>%
  ungroup()

WMH_LowGroup <- WMH_LowGroup %>%
  mutate(Fold = sample(rep(1:n_folds, length.out = n()), n())) %>%
  ungroup()

# Loop through folds
for (fold in 1:n_folds) {
  cat("Processing Fold:", fold, "\n")
  
  # Split into subgroups
  WMH_Low <- WMH_LowGroup %>% filter(WMH_LowGroup$Fold == fold)
  WMH_High <- WMH_HighGroup %>% filter(WMH_HighGroup$Fold == fold)
  
  # Compute adjacency matrices
  Node_WMH_Low <- GenerateNode(WMH_Low)
  Node_WMH_High <- GenerateNode(WMH_High)
  AdjMat_WMH_Low <- Compute_AdjacencyMatrix(Node_WMH_Low)
  AdjMat_WMH_High <- Compute_AdjacencyMatrix(Node_WMH_High)
  
  # Compute thresholded % remaining
  WMH_Low_thr <- extract_thresholded_SCN_data(AdjMat_WMH_Low, thresh_range)
  WMH_High_thr <- extract_thresholded_SCN_data(AdjMat_WMH_High, thresh_range)
  
  # Format data for step function computation and Cox models
  WMH_Low_km <- format_km_data(WMH_Low_thr, paste0("WMH_Low - Fold ", fold))
  WMH_High_km <- format_km_data(WMH_High_thr, paste0("WMH_High - Fold ", fold))
  
  # Store threhold data for this fold
  fold_survival_data <- rbind(WMH_Low_km, WMH_High_km)
  all_survival_data <- rbind(all_survival_data, fold_survival_data)
  
  # (Optional check) Run log-rank test for this fold
  if (length(unique(fold_survival_data$group)) == 2) {
    logrank_test <- survdiff(Surv(time, status) ~ group, data = fold_survival_data)
    chi_sq_value <- logrank_test$chisq
    p_value <- 1 - pchisq(chi_sq_value, df = 1)
    
    # Store results per fold
    logrank_results <- rbind(logrank_results, 
                             data.frame(Fold = fold, Chi_Square = chi_sq_value, P_Value = p_value))
  }
  
  # (Exploratory) Repeat for negative correlations
  NegAdjMat_WMH_Low <- AdjMat_WMH_Low * -1
  NegAdjMat_WMH_High <- AdjMat_WMH_High * -1
  
  NegWMH_Low_survival <- extract_thresholded_SCN_data(NegAdjMat_WMH_Low, thresh_range)
  NegWMH_High_survival <- extract_thresholded_SCN_data(NegAdjMat_WMH_High, thresh_range)
  
  NegWMH_Low_km <- format_km_data(NegWMH_Low_survival, paste0("WMH_Low, Negative - Fold ", fold))
  NegWMH_High_km <- format_km_data(NegWMH_High_survival, paste0("WMH_High, Negative - Fold ", fold))
  
  fold_neg_survival_data <- rbind(NegWMH_Low_km, NegWMH_High_km)
  all_neg_survival_data <- rbind(all_neg_survival_data, fold_neg_survival_data)
  
  # (Exploratory) Run log-rank test for negative correlations per fold
  if (length(unique(fold_neg_survival_data$group)) == 2) {
    neg_logrank_test <- survdiff(Surv(time, status) ~ group, data = fold_neg_survival_data)
    neg_chi_sq_value <- neg_logrank_test$chisq
    neg_p_value <- 1 - pchisq(neg_chi_sq_value, df = 1)
    
    # Store results per fold
    neg_logrank_results <- rbind(neg_logrank_results, 
                                 data.frame(Fold = fold, Chi_Square = neg_chi_sq_value, P_Value = neg_p_value))
  }
}


#------------------------------------------------------------------------------
# 1.5 Apply Cox proportional hazards model to quantify SCN connectivity differences

cox_survival_data <- all_survival_data %>%
  tidyr::separate(group, into = c("Group", "Fold"), sep = " - Fold ") %>%
  mutate(Fold = as.integer(Fold))  # Optional: make Fold a numeric column

cox_survival_data$Group <- relevel(factor(cox_survival_data$Group), ref = "WMH_Low")

cox_model <- coxph(Surv(time, status) ~ Group + strata(Fold), 
                   data = cox_survival_data)

# Show the Cox model results
summary(cox_model)

# plots the scaled Schoenfeld residuals against time for each predictor.
plot(cox.zph(cox_model))
abline(h = 0, lty = 5)


#------------------------------------------------------------------------------
# 1.6 Visualization
# Compute KM fits per fold

all_survival_data <- all_survival_data %>%
  tidyr::separate(group, into = c("WMH_group", "fold"), sep = " - ", remove = FALSE)

km_list <- all_survival_data %>%
  group_split(WMH_group, fold) %>%
  purrr::map(~ {
    survfit(Surv(time, status) ~ 1, data = .x) %>%
      broom::tidy() %>%
      mutate(WMH_group = unique(.x$WMH_group),
             fold = unique(.x$fold))
  })

head(km_list)

# Combine all into one dataframe
km_df <- dplyr::bind_rows(km_list)

# Calculate the 95% Confidence interval via 1.96*SE
km_summary <- km_df %>%
  group_by(WMH_group, time) %>%
  summarise(mean_surv = mean(estimate, na.rm = TRUE),
            sd_surv = sd(estimate, na.rm = TRUE),
            n = n(),
            se_surv = sd_surv / sqrt(n),
            lower = mean_surv - 1.96 * se_surv,
            upper = mean_surv + 1.96 * se_surv,
            .groups = "drop")

# Create step function plots using geom_step function
ggplot(km_summary, aes(x = time, y = mean_surv,
                       color = WMH_group, fill = WMH_group)) +
  geom_step(size = 1.1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, color = NA) +
  labs(x = "Threshold on correlation", y = "Percent of connections remained",
       title = "Group-level Network Survival Curves with 95% CI") +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("WMH_High" = "magenta", "WMH_Low" = "green4")) +
  scale_fill_manual(values = c("WMH_High" = "magenta", "WMH_Low" = "#4DAF4A")) +
  theme(legend.position = "top",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 12),
        axis.title   = element_text(size = 15),
        axis.text    = element_text(size = 13),
        plot.title   = element_text(size = 17, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.3, color = "gray90"))


###############################################################################
# Section 2: To identify candidates of backbone connections within SCNs and 
#            compare the backbone network organization between WMH-high versus
#            WMH-low groups. Finally generate a circular ribbon plot to visualize
#            the change of backbone grey matter network organization.
#
# Summary of methods please refer to the following manuscript figures:
#   1. Supplementary Figure 7
#   2. Method Subheading 5
#   3. Figure 2
###############################################################################
#------------------------------------------------------------------------------
# Function definitions:
# 1. Function to convert correlation matrix to circlize-ready dataframe
corMat_to_df <- function(cor_mat) {
  # Checklist before running the code:
  # Matrix must be square
  stopifnot(nrow(cor_mat) == ncol(cor_mat)) 
  # Matrix must be symmetric with names
  stopifnot(all(rownames(cor_mat) == colnames(cor_mat)))  
  
  rn <- rownames(cor_mat)
  # Grab upper triangle only (no duplicates, no self-loops)
  idx <- which(upper.tri(cor_mat), arr.ind = TRUE)
  
  # Creating the dataframe ready to plug in the circlize library
  df <- data.frame(
    Network1 = rn[idx[,1]],
    Network2 = rn[idx[,2]],
    cor      = cor_mat[idx]
  )
  df
}

# 2.  To create a dataframe of normalized volumetric data. The data frame will
#      be used to generate adjacency matrix for SCN analysis. 
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

# 3. To compute adjacency matrix
Compute_AdjacencyMatrix <- function(data) {
  cor_matrix <- cor(data)                             # Calculate the correlation matrix
  adjacency_matrix <- cor_matrix                      # Create an adjacency matrix
  diag(adjacency_matrix) <- 0                         # Set diagonal to 0 (no self-loops)
  return(adjacency_matrix)
}

# 4. Function to convert the igraph outputs to match with circlize input
to_chord_df <- function(edges,
                        sep = "-",
                        name_map = NULL,      # e.g., c(FOpe = "IFPO")
                        sort_pair = FALSE,    # TRUE if edges are undirected and you want a consistent order
                        drop_duplicates = FALSE) {
  
  sp <- strsplit(edges, split = sep, fixed = TRUE)
  mat <- do.call(rbind, lapply(sp, function(x) {
    x <- trimws(x)
    # keep only first two tokens if there are extras
    if (length(x) < 2) x <- c(x, NA_character_)
    x[1:2]
  }))
  df <- data.frame(Network1 = mat[,1], Network2 = mat[,2], stringsAsFactors = FALSE)
  
  # Standardize names via mapping
  if (!is.null(name_map)) {
    map <- name_map
    df$Network1[df$Network1 %in% names(map)] <- map[df$Network1[df$Network1 %in% names(map)]]
    df$Network2[df$Network2 %in% names(map)] <- map[df$Network2[df$Network2 %in% names(map)]]
  }
  
  # Sort each pair alphabetically for undirected edges
  if (isTRUE(sort_pair)) {
    ord <- t(apply(df, 1, function(r) sort(r)))
    df$Network1 <- ord[,1]
    df$Network2 <- ord[,2]
  }
  
  # Drop duplicate edges (e.g. treating A-B same as B-A if sort_pair=TRUE)
  # If your input graph is directed, i.e. A->B is different than B->A, please
  # Comment out this section!
  if (isTRUE(drop_duplicates)) {
    key <- apply(df, 1, function(r) paste(r, collapse = "||"))
    df <- df[!duplicated(key), , drop = FALSE]
    rownames(df) <- NULL
  }
  
  df
}


# 5. Function to mark vulnerable rich club edges in red color
# note: Sometimes we named backbone "rich edges" but we decided to discard this 
#       name when preparing the manuscript.
mark_richclub_edges <- function(df,
                                missing_edges,
                                preserved_edges,
                                col_vul = "#E31A1C",
                                col_all = "grey80",
                                a_vul = 1.0,  a_all = 0.7, a_other = 0.01,
                                lwd_vul = 20, lwd_all = 20,   lwd_other = 0.5) {
  
  # Check that the main edge list contains the required node columns.
  stopifnot(all(c("Network1","Network2") %in% names(df)))
  names(missing_edges)     <- c("Network1","Network2")
  names(preserved_edges) <- c("Network1","Network2")
  
  clean <- function(x) trimws(as.character(x))
  for (nm in c("Network1","Network2")) {
    df[[nm]]             <- clean(df[[nm]])
    missing_edges[[nm]]     <- clean(missing_edges[[nm]])
    preserved_edges[[nm]] <- clean(preserved_edges[[nm]])
  }
  
  # Create an undirected edge key
  # For undirected graphs only!
  key <- function(a,b) paste(pmin(a,b), pmax(a,b), sep="||")
  df$key          <- key(df$Network1, df$Network2)
  keys_vul        <- unique(key(missing_edges$Network1,     missing_edges$Network2))
  keys_all        <- unique(key(preserved_edges$Network1, preserved_edges$Network2))
  
  # # Mark whether each edge belongs to the preserved backbone edge set or not.
  df$is_all_rich  <- df$key %in% keys_all
  df$is_vul       <- df$key %in% keys_vul
  
  # colors with baked alpha
  df$col <- adjustcolor("grey90", alpha.f = a_other)
  df$col[df$is_all_rich] <- adjustcolor(col_all, alpha.f = a_all)
  df$col[df$is_vul]      <- adjustcolor(col_vul, alpha.f = a_vul)
  
  # absolute widths
  df$lwd <- lwd_other
  df$lwd[df$is_all_rich] <- lwd_all
  df$lwd[df$is_vul]      <- lwd_vul
  
  df
}

#------------------------------------------------------------------------------
# 2.1 Parameter input and abbreviation definition

# Number of replications for concensus matrix
B <- 100        

# choose stability threshold (e.g., must appear in ≥50% of bootstraps)
stable_cut <- 0.5

# Folds of data per group
n_folds <- 8          

# Percentile threshold for an edge to be considered as backbone edge
Backbone_threshold <- 0.8

# Example default setting interpretations: 
# Backbone threshold=0.8, stable cut=0.5, B=100, n_fold=8
# To be considered as a part of the backbone network of grey matter SCN, that 
# edge in the SCN graph needs to be ranked as top 20% in terms of its weights 
# in all 8 folds of data, and show up at least 50% times among 100 repeats of 
# random fold generation.

# Create a list of abbreviations
abbreviation_table <- c(
  "aCG" = "Anterior Cingulate Gyrus",
  "Amyg" = "Amygdala",
  "COC" = "Central Opercular Cortex",
  "MFG" = "Middle Frontal Gyrus",
  "NAcc" = "Nucleus Accumbens",
  "PaCG" = "Paracingulate Gyrus",
  "Pall" = "Pallidum",
  "POC" = "Parietal operculum",
  "SFG" = "Superior Frontal Gyrus",
  "Scal" = "Subcallosal Cortex",
  "Thal" = "Thalamus",
  "pCG" = "Posterior Cingulate Gyrus",
  "aPHG" = "Anterior Parahippocampal Gyrus",
  "pITG" = "Posterior Inferior Temporal Gyrus",
  "pMTG" = "Posterior Middle Temporal Gyrus",
  "toMTG" = "Temporooccipital Middle Temporal Gyrus",
  "aITG" = "Anterior Inferior Temporal Gyrus",
  "PP" = "Planum Polare",
  "FMC" = "Middle Frontal Cortex",
  "Hipp" = "Hippocampus",
  "Puta" = "Putamen",
  "IFPT" = "Inferior Frontal Triangularis"
)

group_lookup <- c(
  # Frontal
  "PaCG"="Frontal","OFC"="Frontal","FP"="Frontal","FOpe"="Frontal",
  "MFG"="Frontal","SFG"="Frontal","FMC"="Frontal","IFPT"="Frontal","IFPO"="Frontal",
  # Parietal
  "Prec"="Parietal","AG"="Parietal","aSMG"="Parietal",
  "pSMG"="Parietal","POC"="Parietal","COC"="Parietal",
  # Temporal
  "aPHG"="Temporal","pPHG"="Temporal",
  "aITG"="Temporal","pITG"="Temporal","pMTG"="Temporal",
  "toMTG"="Temporal","aMTG"="Temporal","aSTG"="Temporal",
  "pSTG"="Temporal","toITG"="Temporal","PP"="Temporal",
  # Occipital
  
  # Subcortex
  "Thal"="Subcortex","Caud"="Subcortex","Puta"="Subcortex","Pall"="Subcortex",
  "NAcc"="Subcortex","Hipp"="Subcortex","Amyg"="Subcortex",
  # Limbic/Other (adjust if needed)
  "Scal"="Other","Insu"="Other","aCG"="Other", "pCG"="Other"
)

#------------------------------------------------------------------------------
# 2.2 Detect backbone edges
# Creating result data frames
boot_shared_list    <- vector("list", B)  # shared rich edges per iteration
boot_vulnerable_list <- vector("list", B) # vulnerable edges per iteration

# Loop through B iterations to search for backbone edges
for (b in 1:B) {
  cat("Iteration", b, "of", B, "\n")
  
  # Fold assignment without replacement
  older_boot <- older %>%
    group_by(WMH_status) %>%
    group_modify(~ {
      group_size <- nrow(.x)
      .x[sample(seq_len(group_size), size = group_size, replace = FALSE), ]
    }) %>%
    ungroup()
  
  # Randomly assign folds within the sample:
  older_boot <- older_boot %>%
    mutate(Fold = sample(rep(1:n_folds, length.out = n()), n())) %>%
    ungroup()
  
  # Tracking top edges per fold
  edge_appearance_H <- vector("list", n_folds)
  edge_appearance_L <- vector("list", n_folds)
  
  for (fold in 1:n_folds) {
    # split by WMH group + fold
    High <- older_boot %>% filter(WMH_status == "High", Fold == fold)
    Low  <- older_boot %>% filter(WMH_status == "Low",  Fold == fold)
    
    # skip empty folds just in case
    if (nrow(High) < 3 || nrow(Low) < 3) next
    
    # adjacency matrices
    Adj_H <- Compute_AdjacencyMatrix(GenerateNode(High))
    Adj_L <- Compute_AdjacencyMatrix(GenerateNode(Low))
    
    # rank edges within each fold
    df_H <- get_edge_df(Adj_H, "WMH_High") %>%
      mutate(Rank = rank(Weight) / n())
    df_L <- get_edge_df(Adj_L, "WMH_Low") %>%
      mutate(Rank = rank(Weight) / n())
    
    # top edges (top 20%)
    top_H_edges <- df_H$Edge[df_H$Rank >= Backbone_threshold]
    top_L_edges <- df_L$Edge[df_L$Rank >= Backbone_threshold]
    
    edge_appearance_H[[fold]] <- top_H_edges
    edge_appearance_L[[fold]] <- top_L_edges
  }
  
  # Count edges that being marked as backbone edges in all folds of data
  all_H <- unlist(edge_appearance_H)
  all_L <- unlist(edge_appearance_L)
  
  if (length(all_H) == 0 || length(all_L) == 0) {
    boot_shared_list[[b]]     <- character(0)
    boot_vulnerable_list[[b]] <- character(0)
    next
  }
  
  tab_H  <- table(all_H)
  tab_L  <- table(all_L)
  
  rich_H_b <- names(tab_H)[tab_H == n_folds]
  rich_L_b <- names(tab_L)[tab_L == n_folds]
  
  # shared & vulnerable for run
  shared_rich_b <- intersect(rich_H_b, rich_L_b)
  Vulnerable_b  <- setdiff(rich_L_b, rich_H_b)   # present in Low only
  
  boot_shared_list[[b]]     <- shared_rich_b
  boot_vulnerable_list[[b]] <- Vulnerable_b
}


# Report the observed backbone edges in the SCNs:
# Aggregate over B iterations: keep *stable* edges
# how often each edge is shared/vulnerable across iterations
all_shared    <- unlist(boot_shared_list)
all_vulnerable <- unlist(boot_vulnerable_list)

# Calculate the frequency which edges are "nominated" as backbone connections.
# Shared edges are backbone connections observed in both WMH-low and WMH-high group
# Vulnerable edges are backbone connections only observed in WMH-low group.
freq_shared    <- table(all_shared)    / B
freq_vulnerable <- table(all_vulnerable) / B

print(freq_shared)
print(freq_vulnerable)

shared_rich <- names(freq_shared)[freq_shared >= stable_cut]
Vulnerable  <- names(freq_vulnerable)[freq_vulnerable >= stable_cut]

cat("Number of shared backbone edges:", length(shared_rich), "\n")
cat("Number of WMH-vulnerable backbone edges:", length(Vulnerable),  "\n")


#------------------------------------------------------------------------------
# 2.3 Visualize the results using igraph and ribbon plot

# Collectively combine the required data for plotting
all_edges <- unique(c(shared_rich, Vulnerable))
edge_list <- do.call(rbind, strsplit(as.character(all_edges), "-"))
g_rich_all <- graph_from_edgelist(as.matrix(edge_list), directed = FALSE)

# Create an igraph network object from an edge list
edge_names_graph <- apply(ends(g_rich_all, E(g_rich_all)), 1, function(x) paste(sort(x), collapse = "-"))
shared_set <- sort(unique(as.character(shared_rich)))
Vulnerable_set <- sort(unique(as.character(Vulnerable)))

# Initial glance on SCN graphs. Vulnerable edges are marked red in the diagram
edge_colors <- ifelse(edge_names_graph %in% Vulnerable_set, "red",
                      ifelse(edge_names_graph %in% shared_set, "grey", "grey"))

fixed_layout_wmh_diff <- layout_with_fr(g_rich_all)
plot(g_rich_all,
     vertex.label = V(g_rich_all)$name,
     edge.color = edge_colors,
     edge.width = 5,
     vertex.size = 15,
     vertex.label.cex = 0.8,
     vertex.label.color = "black",
     layout = fixed_layout_wmh_diff,
     main = "Rich Club Network\nGrey = Shared, Red = Lost Connection, Green = Gained Connection")

# Initializing the dataframe and parameter input to circlize library

df_chord <- corMat_to_df(Compute_AdjacencyMatrix(GenerateNode(Low)))

# Circlize appearance settings:
tau   <- 0.1   # threshold where ramp gets steep
gamma <- 2     # >1 makes post-tau ramp steeper; try 1.4–2.0

# Apply the customized appearance settings to circlize plotting function
df <- df_chord %>%
  mutate(
    # only visualize positive correlations
    val_pos = pmax(cor, 0),
    # piecewise squash: 0..tau -> ~0; tau..1 -> rescaled to 0..1
    v_lin = pmin(pmax((val_pos - tau) / (1 - tau), 0), 1),
    # optional gamma to steepen further
    v = v_lin ^ gamma,
    # Color: greyscale with stronger edges near black only when v is high
    col = colorRamp2(c(0, 1), c("white", "black"))(v),
    # Alpha: mostly faint below tau, opaque after
    alpha = rescale(v, to = c(0.05, 1)),
    # Line width: thin below tau, thicker after - default c(0.2, 10)
    lwd = rescale(v, to = c(0.2, 10))
  )

# Create sector labels based on the brain region definitions entered in line 601
all_nodes <- sort(unique(c(df$Network1, df$Network2)))

roi_groups <- tibble(
  ROI   = all_nodes,
  Group = unname(group_lookup[all_nodes])
) %>%
  mutate(
    Group = ifelse(is.na(Group), "NA", Group),
    # coalesce subgroups -> main lobe groups
    Group = case_when(
      str_detect(Group, "^Frontal")   ~ "Frontal",
      str_detect(Group, "^Temporal")  ~ "Temporal",
      str_detect(Group, "^Parietal")  ~ "Parietal",
      str_detect(Group, "^Subcortex") ~ "Subcortex",
      str_detect(Group, "^Other") ~ "Other",
      TRUE ~ Group
    ),
    Group = factor(Group, levels = c("Frontal","Parietal","Temporal","Subcortex","Other"))
  ) %>%
  arrange(Group, ROI)
sector_order <- roi_groups$ROI

# Transform the igraph ouput into the dataframe useable by chord function: 
Vulnerable_edges <- to_chord_df(
  Vulnerable,
  sort_pair = FALSE,     # set TRUE if you want undirected canonical ordering
  drop_duplicates = TRUE
)

backbone_edges <- to_chord_df(
  shared_rich,
  sort_pair = FALSE,     # set TRUE if you want undirected canonical ordering
  drop_duplicates = TRUE
)

# Inspect the vulnerable vs. backbone edges
print(Vulnerable_edges)
print(backbone_edges)


# Apply to the plotting dataframe
df_marked <- mark_richclub_edges(
  df,
  missing_edges   = Vulnerable_edges,
  preserved_edges = backbone_edges
)


# Rank the pearson correlations
df_marked <- df_marked %>%
  mutate(
    rank_tmp = ifelse(is_all_rich, cor, NA_real_),     # only keep cor for rich
    rank_rich = dense_rank(rank_tmp),                 # lowest cor = 1, 2, ...
    rank_final = ifelse(is_all_rich, rank_rich + 1, 1) # shift by +1, non-rich = 1
  )

df_marked$rank_final <- df_marked$rank_final / (sum(df_marked$is_all_rich) + 2)

# Tabular output of pre-plot backbone edge information.
head(df_marked)
write.csv(df_marked, "BackboneAnalyis_Results.csv")

# Initiate the circos plotting parameters: 
circos.clear()
small_gap <- 1
big_gap   <- 4
runs <- rle(as.character(roi_groups$Group))
ends <- cumsum(runs$lengths)
gap.after <- rep(small_gap, length(sector_order))
gap.after[ends] <- big_gap
group_ht <- 0.06   # middle band: group label track
roi_ht   <- 0.01   # outer band: ROI names (default 0.06)
sep_ht   <- 0.02   # inner thin band: grey separator per ROI (default 0.02)

circos.par(
  start.degree = 90,
  gap.after    = gap.after,
  track.margin = c(0.001, 0.001), # Default 0.001 for both
  cell.padding = c(0,0,0,0)
)

group_label_cols <- c(
  "Frontal"   = "#1f77b4",
  "Parietal"  = "#2ca02c",
  "Temporal"  = "#ff7f0e",
  "Subcortex" = "#8c564b",
  "Other"     = "#9467bd"
)

# build chords and pre-allocate 3 tracks
#   Track 1 = separators (inner), 
#   Track 2 = group labels (middle), 
#   Track 3 = ROI names (outer)
chordDiagram(
  x = df_marked[, c("Network1","Network2","val_pos")],
  order              = sector_order,
  #grid.col           = grid.col,
  col                = df_marked$col,
  link.lwd           = df_marked$rank_final,
  link.sort          = TRUE,
  link.largest.ontop = TRUE,
  transparency       = 0,
  annotationTrack    = NULL,
  preAllocateTracks  = list(
    list(track.height = sep_ht),
    list(track.height = group_ht),
    list(track.height = roi_ht)
  )
)

# Customize the Track 1 (inner) - Showing grey separator per ROI
circos.trackPlotRegion(track.index = 1, bg.border = NA, panel.fun = function(x, y) {
  xl <- CELL_META$xlim
  yl <- CELL_META$ylim
  # fill each sector's strip
  circos.rect(xl[1], yl[1], xl[2], yl[2], col = "grey93", border = NA)
  # subtle border to delineate ROIs
  circos.rect(xl[1], yl[1], xl[2], yl[2], col = NA, border = "grey85", lwd = 1)
})


# Customize the outer track: group boxes + text
# Note: On the plot window the texts show huge size, but as we save the figure
#       the font size scales based on the figure size
table(roi_groups$Group)
grp_bounds <- split(roi_groups$ROI, roi_groups$Group)
for (g in names(grp_bounds)) {
  highlight.sector(
    grp_bounds[[g]],
    track.index = 2,
    col    = NA,
    border = "grey70", lwd = 1,
    text   = g, niceFacing = TRUE, cex = 4.5,
    text.col = group_label_cols[g]
  )
  
# Add ROI names on outermost band
label_outer_mm <- 8
circos.trackPlotRegion(track.index = 3, bg.border = NA, panel.fun = function(x, y) {
    circos.text(
      x = CELL_META$xcenter,
      y = CELL_META$ylim[2] + mm_y(label_outer_mm),
      labels = CELL_META$sector.index,
      facing = "clockwise", niceFacing = TRUE,
      adj = c(0, 0.5), cex = 3.5
    )
  })
}

# The font is set for 3500x3500+ pixel figure save
png("Circular_RichClub.png",
    width = 3500, height = 3500, res = 400)
