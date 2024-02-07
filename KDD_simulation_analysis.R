library(ggplot2)
library(reshape2)
library(rBeta2009)

###################################################################
####################### Part I: Preparations ######################
###################################################################

# params
N_g = 1000
n = 1000
n_iter = 100
index = 1:n
trunc_level = 10

# true outcome model
outcome_Y <- function(M, X){
  lambda = 0.5*sqrt(M) + 2*X + ifelse(X > 0.5, 0.5*X, 0) +  0.5*sqrt(M)*X
  return(rexp(1, rate = 1/lambda))
}

sim_x = rbeta(1000, 0.5, 0.4)
estimate_Y = apply(as.matrix(sim_x), 1, function(x) outcome_Y(0, x))
plot(sim_x, estimate_Y)

# saving
tT_counts_save = cT_counts_save = matrix(NA, ncol = trunc_level, nrow = n_iter)
cC_counts_save = rep(NA, n_iter)

true_tau_M_save = matrix(NA, ncol = trunc_level, nrow = n_iter)
bias_baseline = matrix(NA, ncol = trunc_level, nrow = n_iter)
bias_xgb_class = matrix(NA, ncol = trunc_level, nrow = n_iter)
bias_xgb_class_H = matrix(NA, ncol = trunc_level, nrow = n_iter)
bias_glm = matrix(NA, ncol = trunc_level, nrow = n_iter)
bias_glm_H = matrix(NA, ncol = trunc_level, nrow = n_iter)
bias_treat_minus_controlcontrol = matrix(NA, ncol = trunc_level, nrow = n_iter)
bias_treat_only = matrix(NA, ncol = trunc_level, nrow = n_iter)
baseline_save = treat_minus_controlcontrol_save = treat_only_save = matrix(NA, ncol = trunc_level, nrow = n_iter)
xgb_class_H_save = glm_H_save = matrix(NA, ncol = trunc_level, nrow = n_iter)



# Settings of interference 
# set I
prob_treat_per_game = c(0.4, 0.1, 0.1, 0.1, 0.1, 0.2) 
N_g = 2000
# set II
prob_treat_per_game = c(0.06, 0.02, 0.19, 0.23, 0.34, 0.16) 
N_g = 1000
# set III
prob_treat_per_game = c(0.20, 0.34, 0.07, 0.16, 0.06, 0.17) 
N_g = 1000





###################################################################
################ Part II: Fit Model and Estimation ################
###################################################################
for (iter in 1:n_iter){
  ##########################   ##########################   ########################## 
  ##################### Generate True Data Frame and Obtain True Tau M ###############
  ##########################   ##########################   ########################## 
  print(paste("The current iter is:", iter))
  rm(df)
  
  # observed data frame df
  ## pre-experiment: pre_X and pre_exp Y0 for all players
  pre_X = rbeta(n, 0.5, 0.5)
  preExp_Y0 = rep(NA, n)
  for (nn in 1:n){
    preExp_Y0[nn] = outcome_Y(0, pre_X[nn])
  }
  # use preExp_Y0 as unbiased estimator of mu hat（base effect）
  
  ## experiment: Z, M, Y0 and Yt
  Z = rbinom(n, 1, 0.5) # randomized Z
  
  ### play games
  # Create adjusted probabilities of X
  # Normalize x
  normalized_x <- pre_X / sum(pre_X)
  # Here we assume a base probability of 0.8 (to complement the 0.2 correlation with normalized_x)
  adjusted_prob_X <- 0.8 * (1 / length(pre_X)) + 0.2 * normalized_x^2
  # Ensure the probabilities sum to 1
  adjusted_prob_X <- adjusted_prob_X / sum(adjusted_prob_X)
  
  X_Z1 = normalized_x[which(Z == 1)]
  X_Z0 = pre_X[which(Z == 0)]
  adjusted_prob_X_Z1 = 0.8 * (1 / length(X_Z1)) + 0.2 * X_Z1^2
  adjusted_prob_X_Z0 = adjusted_prob_X[which(Z == 0 & pre_X >= 0.2)]
  # adjusted_prob_X_Z0 = 0.8 * (1 / length(X_Z0[which(X_Z0 >= 0.2)])) + 0.2 * X_Z0[which(X_Z0 >= 0.2)]^2
  nonadjusted_prob_X_Z0 = pre_X[which(Z == 0 & pre_X < 0.2)]
  
  # number of games
  M_cT = rep(0, n)
  M_cC = rep(0, n)
  M_tT = rep(0, n)
  
  for (j in 1:N_g) {
    num_treat_per_game = sample(0:5, 1, replace = FALSE, prob = prob_treat_per_game)
    treat_index = sample(index[which(Z == 1)], num_treat_per_game, prob = adjusted_prob_X_Z1 / sum(adjusted_prob_X_Z1))
    if (num_treat_per_game == 0){
      control_index = sample(index[which(Z == 0 & pre_X < 0.2)], 5, prob = nonadjusted_prob_X_Z0 / sum(nonadjusted_prob_X_Z0))
    }else{
      control_index = sample(index[which(Z == 0 & pre_X >= 0.2)], 5 - num_treat_per_game, prob =  adjusted_prob_X_Z0 / sum(adjusted_prob_X_Z0))
    }
    player_index = c(treat_index, control_index)
    
    for (i in 1:5){
      idx = player_index[i] # current player index i in jth game
      Z_ij = Z[player_index][i] # current player original treatment assignment
      phi = (sum(Z[player_index]) - Z_ij) / (5-1) # the proportion of the team member in treatment
      if (Z_ij == 1){
        M_tT[idx] = M_tT[idx] + 1
      }
      if (Z_ij == 0){
        if (phi > 0) {
          M_cT[idx] = M_cT[idx] + 1
        } else {
          M_cC[idx] = M_cC[idx] + 1
        }
      }
    }
  }
  
  # combine df as observed data frame
  df = cbind(index, pre_X, preExp_Y0, Z, M_tT, M_cT, M_cC)
  df = as.data.frame(df)
  
  # calculate Yt for treat and control_mixed players
  Yt = rep(NA, n)
  
  # the value of Ms for the treatment games
  Ms = 1:max(as.numeric(names(table(df$M_tT[df$Z == 1]))))
  true_tau_M = rep(NA, length(Ms))
  
  for (m in 1:length(Ms)){
    M = Ms[m]
    YM = rep(0, n)
    for (j in 1:n){
      YM[j] = outcome_Y(M, df$pre_X[j])
      
      if (df$Z[j] == 1 & df$M_tT[j] == M){
        Yt[j] = YM[j]
      }else if (df$Z[j] == 0 & df$M_cT[j] == M){
        Yt[j] = YM[j]
      }
    }
    # calculate true tau_M as bar(Y(M)) - bar(Y(0))
    mean_YM = mean(YM)
    mean_Y0 = mean(preExp_Y0) 
    true_tau_M[m] = mean_YM - mean_Y0
  }
  # save yt to the df
  df$Yt = Yt
  
  # check if any treat player plays no games at all
  df[which(df$M_tT == 0 & df$Z == 1),]
  # check the control-control players
  df[which(df$M_cT == 0 & df$Z == 0),]
  
  tT_counts = table(df$M_tT[which(df$Z == 1)])
  tT_counts = as.data.frame(tT_counts)
  cT_counts = table(df$M_cT[which(df$Z == 0 & df$M_cT > 0)])
  cT_counts = as.data.frame(cT_counts)
  
  for (m in 1:trunc_level){
    if (m < trunc_level){
      tT_counts_save[iter, m] = length(df$M_tT[which(df$Z == 1 & df$M_tT == m)])
      cT_counts_save[iter, m] = length(df$M_cT[which(df$Z == 0 & df$M_cT == m)])
    }else{
      tT_counts_save[iter, m] = length(df$M_tT[which(df$Z == 1 & df$M_tT >= trunc_level)])
      cT_counts_save[iter, m] = length(df$M_cT[which(df$Z == 0 & df$M_cT >= trunc_level)])
    }
  }
  cC_counts_save[iter] = dim(df[which(df$Z == 0 & df$M_cT == 0), ])[1]
  
  
  
  ##########################   ##########################   ########################## 
  ########################## Calculate Estimated Tau M #############################
  ##########################   ##########################   ########################## 
  library(stats)
  library(MASS)
  library(xgboost)
  library(rBeta2009)
  
  # mu hat estimation
  df$predicted_Y0 = preExp_Y0
  
  # generalized propensity score
  df$M = df$M_tT + df$M_cT
  
  
  ########### XgBoost as Classification ###########
  # truncate at trunc_level
  df$M_truncated = ifelse(df$M >= trunc_level, trunc_level, df$M)
  table(df$M_truncated)
  
  # Prepare data for XGBoost
  data_matrix <- xgb.DMatrix(as.matrix(cbind(df$pre_X, df$Z)), label = as.numeric(as.factor(df$M_truncated)) - 1)
  num_class <- length(unique(df$M_truncated))
  # Set parameters for XGBoost
  params <- list(
    objective = "multi:softprob",
    num_class = num_class,
    eta = 0.5,
    max_depth = 6
  )
  
  # Train the model
  xgb_class_model <- xgboost(data = data_matrix, params = params, nrounds = 500)
  
  # Predict probabilities
  pred_prob <- predict(xgb_class_model, data_matrix)
  
  # Reshape probabilities to a matrix (columns as classes)
  prob_matrix <- matrix(pred_prob, nrow = length(pred_prob)/num_class, byrow = TRUE)
  
  # Extract probabilities of the observed classes
  gps_xgb_class <- prob_matrix[cbind(1:nrow(prob_matrix), as.numeric(as.factor(df$M_truncated)))]
  
  df$gps_xgb_class = gps_xgb_class
  
  
  
  
  # tau M estimation
  df_treat_and_control_mixed = df[which(df$Z != 0 | df$M_cT > 0), ]
  df_treat = df[which(df$Z == 1), ]
  df_control_mixed = df[which(df$Z == 0 & df$M_cT > 0), ]
  df_control_control = df[which(df$Z == 0 & df$M_cT == 0), ]
  # truncate
  df_treat$Mt_t_truncated = ifelse(df_treat$M_tT >= trunc_level, trunc_level, df_treat$M_tT)
  df_control_mixed$Mc_t_truncated = ifelse(df_control_mixed$M_cT >= trunc_level, trunc_level, df_control_mixed$M_cT)
  
  M_table = table(c(df_treat$Mt_t_truncated, df_control_mixed$Mc_t_truncated))
  Ms = as.numeric(names(M_table))
  
  sim_predicted_tau_M_xgb_class_H = rep(NA, length(Ms))
  for (i in 1:length(Ms)){
    m = Ms[i]
    sim_predicted_tau_M_xgb_class_H[i] = (sum(df_treat$Yt[which(df_treat$Mt_t_truncated == m)] / df_treat$gps_xgb_class[which(df_treat$Mt_t_truncated == m)]) + sum(df_control_mixed$Yt[which(df_control_mixed$Mc_t_truncated == m)] / df_control_mixed$gps_xgb_class[which(df_control_mixed$Mc_t_truncated == m)])) /  
      (sum(ifelse(df_treat$Mt_t_truncated == m, 1/df_treat$gps_xgb_class, 0)) + sum(ifelse(df_control_mixed$Mc_t_truncated == m, 1/df_control_mixed$gps_xgb_class, 0))) -
      (sum(df_treat$predicted_Y0) + sum(df_control_mixed$predicted_Y0) + sum(df_control_control$predicted_Y0)) / (dim(df)[1])
  }
  
  
  
  ################### Baseline Tau ############################
  tau_hat_1 = rep(NA, length(unique(M_table)))
  tau_baseline =  rep(NA, length(unique(M_table)))
  tau_M_hat_treat = rep(NA, length(Ms))
  for (i in 1:length(Ms)){
    M = Ms[i]
    tau_baseline[i] = mean(df_treat$Yt[which(df_treat$M_tT == M)]) - mean(c(df_control_mixed$Yt, df_control_control$preExp_Y0))
    tau_hat_1[i] = mean(df_treat$Yt[which(df_treat$M_tT == M)]) - mean(df_control_control$preExp_Y0) 
    tau_M_hat_treat[i] = mean(df_treat$Yt[which(df_treat$M_tT == M)]) - mean(df_treat$predicted_Y0)
  }
  
  ###### obtain bias of tau M ###### 
  bias_xgb_class_H[iter, Ms[which(Ms >= 1 & Ms <= trunc_level)]] = true_tau_M[which(Ms >= 1 & Ms <= trunc_level)] - sim_predicted_tau_M_xgb_class_H[which(Ms >= 1 & Ms <= trunc_level)]
  bias_baseline[iter, Ms[which(Ms >= 1 & Ms <= trunc_level)]] = true_tau_M[which(Ms >= 1 & Ms <= trunc_level)] - tau_baseline[which(Ms >= 1 & Ms <= trunc_level)]
  bias_treat_minus_controlcontrol[iter, Ms[which(Ms >= 1 & Ms <= trunc_level)]] = true_tau_M[which(Ms >= 1 & Ms <= trunc_level)] - tau_hat_1[which(Ms >= 1 & Ms <= trunc_level)]
  bias_treat_only[iter, Ms[which(Ms >= 1 & Ms <= trunc_level)]] = true_tau_M[which(Ms >= 1 & Ms <= trunc_level)] - tau_M_hat_treat[which(Ms >= 1 & Ms <= trunc_level)]
  
  ###### save true tau M and estimated tau M ######
  true_tau_M_save[iter, Ms[which(Ms >= 1 & Ms <= trunc_level)]] = true_tau_M[which(Ms >= 1 & Ms <= trunc_level)]
  baseline_save[iter, Ms[which(Ms >= 1 & Ms <= trunc_level)]] = tau_baseline[which(Ms >= 1 & Ms <= trunc_level)]
  treat_minus_controlcontrol_save[iter, Ms[which(Ms >= 1 & Ms <= trunc_level)]] = tau_hat_1[which(Ms >= 1 & Ms <= trunc_level)]
  treat_only_save[iter, Ms[which(Ms >= 1 & Ms <= trunc_level)]] = tau_M_hat_treat[which(Ms >= 1 & Ms <= trunc_level)]
  xgb_class_H_save[iter, Ms[which(Ms >= 1 & Ms <= trunc_level)]] = sim_predicted_tau_M_xgb_class_H[which(Ms >= 1 & Ms <= trunc_level)]
}



######################################################################################
########################## Part III: Visualization ###################################
######################################################################################

#########################  Simulation visualization for interference ##############
# Convert the matrices to data frames
cT_counts_df <- data.frame(M = c(1:9, "10+"), t(cT_counts_save))
tT_counts_df <- data.frame(M = c(1:9, "10+"), t(tT_counts_save))
cC_counts_df <- data.frame(M = 0, t(cC_counts_save))
# Add a column 'Group' to distinguish between cT_counts_save and tT_counts_save
cT_counts_df$Group <- 'Control_Mixed'
tT_counts_df$Group <- 'Treat'
cC_counts_df$Group <- 'Control_Control'

# Combine the two data frames
combined_df <- rbind(cT_counts_df, tT_counts_df, cC_counts_df)
combined_df$M <- factor(combined_df$M, levels = c(0, 1:9, "10+"))

# Melt the data frame to long format
combined_df_long <- melt(combined_df, id.vars = c("M", "Group"), variable.name = "Iteration", value.name = "Frequency")
combined_df_long$Iteration = as.factor(combined_df_long$Iteration)
# Create the ggplot
df_descrip_plot = ggplot(combined_df_long, aes(x = M, y = Frequency/n, group = interaction(Iteration,Group), color = Group)) +
  geom_line(alpha = 0.1) +
  geom_point(size = 0.8, alpha = 0.1)+
  theme_minimal() +
  labs(x = "M", y = "Sample Size Proportion", color = "Group") +
  scale_color_manual(values = c("Control_Mixed" = "#2166AC", "Treat" = "#B2182B", "Control_Control" = "burlywood"))
df_descrip_plot

df_descrip_Set_I_update = ggplot(combined_df_long, aes(x = M, y = Frequency/n, group = interaction(Iteration, Group), color = Group)) +
  geom_line(aes(alpha = ifelse(Group %in% c("Control_Mixed", "Treat"), 0.1, 1))) +
  geom_point(aes(alpha = ifelse(Group %in% c("Control_Mixed", "Treat"), 0.1, 1)), size = 0.8) +
  theme_minimal() +
  labs(x = "M", y = "Sample Size Proportion", color = "Group") +
  scale_color_manual(values = c("Control_Mixed" = "#2166AC", "Treat" = "#B2182B", "Control_Control" = "burlywood")) +
  theme(text = element_text(size = 24,face = "bold"),  
        legend.title = element_blank(), 
        legend.text=element_text(size = 25,face = "bold"), 
        axis.title = element_text(size = 24, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"))  +
  theme(
    legend.position = "bottom",  
    legend.box = "horizontal",   
    legend.box.just = "center",  
    legend.margin = margin(6, 6, 6, 6)  # Adjust the margin around the legend
  ) +
  guides(
    alpha = "none", # Disable the alpha guide
    color = guide_legend(override.aes = list(alpha = 1, size = 4, shape = 16)),
    fill = guide_legend(override.aes = list(size = 2))# Override alpha for color legend to be fully opaque
  )




######################### bias boxplot ###########################
df_bias_treat_minus_controlcontrol <- data.frame(bias_treat_minus_controlcontrol)
df_bias_treat_minus_controlcontrol$Method <- 'Naive - w/o Control-Mixed'

df_bias_treat_only <- data.frame(bias_treat_only)
df_bias_treat_only$Method <- 'Proposed - w/o Control-Mixed'

df_bias_xgb_class_H <- data.frame(bias_xgb_class_H)
df_bias_xgb_class_H$Method <- 'Proposed'

df_bias_baseline <- data.frame(bias_baseline)
df_bias_baseline$Method <- 'Naive'

# Combine the two data frames
combined_df <- rbind(df_bias_xgb_class_H, df_bias_baseline, df_bias_treat_minus_controlcontrol, df_bias_treat_only)
colnames(combined_df) = c(c(1:9, "10+"), 'Method')

# Convert 'Method' to a factor with levels in the desired order
combined_df$Method <- factor(combined_df$Method, levels = c("Naive", "Naive - w/o Control-Mixed", "Proposed - w/o Control-Mixed", "Proposed"))

# Melt the combined data frame to long format
bias_df <- melt(combined_df, id.vars = 'Method', variable.name = 'M', value.name = 'Bias')

# Create the boxplot using ggplot2
bias_plot_I = ggplot(bias_df, aes(x = M, y = Bias, fill = Method)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Naive" = "darkgrey", "Naive - w/o Control-Mixed" = "lightgrey", "Proposed - w/o Control-Mixed" = "cadetblue",  "Proposed" = "coral3")) +
  theme_minimal() +
  ylim(-8, 8) + 
  labs(x = "M", y = "Bias")
bias_plot_I


######################### line plots for estimation tau M comparisons ######################### 
df_treat_minus_controlcontrol_save <- data.frame(mean =  colMeans(treat_minus_controlcontrol_save, na.rm = TRUE),
                                                 lower = apply(treat_minus_controlcontrol_save, 2, function(x) quantile(x, 0.025, na.rm = TRUE)),
                                                 upper = apply(treat_minus_controlcontrol_save, 2, function(x) quantile(x, 0.975, na.rm = TRUE)),
                                                 Method = rep('Naive - w/o Control-Mixed', trunc_level),
                                                 truth = colMeans(true_tau_M_save, na.rm = TRUE),
                                                 M = c(1:9, "10+"))


df_treat_only_save <- data.frame(mean =  colMeans(treat_only_save, na.rm = TRUE),
                                 lower = apply(treat_only_save, 2, function(x) quantile(x, 0.025, na.rm = TRUE)),
                                 upper = apply(treat_only_save, 2, function(x) quantile(x, 0.975, na.rm = TRUE)),
                                 Method = rep('Proposed - w/o Control-Mixed', trunc_level),
                                 truth = colMeans(true_tau_M_save, na.rm = TRUE),
                                 M = c(1:9, "10+"))

df_xgb_class_H_save <- data.frame(mean =  colMeans(xgb_class_H_save, na.rm = TRUE),
                                  lower = apply(xgb_class_H_save, 2, function(x) quantile(x, 0.025, na.rm = TRUE)),
                                  upper = apply(xgb_class_H_save, 2, function(x) quantile(x, 0.975, na.rm = TRUE)),
                                  Method = rep('Proposed', trunc_level),
                                  truth = colMeans(true_tau_M_save, na.rm = TRUE),
                                  M = c(1:9, "10+"))

df_baseline_save <- data.frame(mean =  colMeans(baseline_save, na.rm = TRUE),
                               lower = apply(baseline_save, 2, function(x) quantile(x, 0.025, na.rm = TRUE)),
                               upper = apply(baseline_save, 2, function(x) quantile(x, 0.975, na.rm = TRUE)),
                               Method = rep('Naive', trunc_level),
                               truth = colMeans(true_tau_M_save, na.rm = TRUE),
                               M = c(1:9, "10+"))


# Combine the two data frames
combined_df_save <- rbind(df_xgb_class_H_save, df_baseline_save, df_treat_minus_controlcontrol_save, df_treat_only_save)

# Convert 'Method' to a factor with levels in the desired order
combined_df_save$Method <- factor(combined_df_save$Method, levels = c("Truth", "Naive", "Naive - w/o Control-Mixed", "Proposed - w/o Control-Mixed", "Proposed"))
combined_df_save$M <- factor(combined_df_save$M, levels = c(1:9, "10+"))


# Create the boxplot using ggplot2
gg_combined_df_save_I_updated = ggplot(combined_df_save, aes(y=mean, x=M, group=Method)) +
  geom_line(aes(color = Method), size = 1)+
  geom_point(aes(color = Method), size = 3)+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Method), alpha = 0.5) +
  geom_line(aes(x = M, y = truth, color = "Truth"),size = 1, alpha = 0.7) +
  geom_point(aes(x = M, y = truth, color = "Truth"), size = 2, alpha = 0.7) +
  facet_wrap(~Method, nrow = 1)+
  scale_fill_manual(values = c("Truth" = "black", "Naive" = "darkgrey", "Naive - w/o Control-Mixed" = "lightgrey","Proposed - w/o Control-Mixed" = "cadetblue", "Proposed" = "coral3")) +
  scale_color_manual(values = c("Truth" = "black", "Naive" = "darkgrey", "Naive - w/o Control-Mixed" = "lightgrey","Proposed - w/o Control-Mixed" = "cadetblue", "Proposed" = "coral3")) +
  theme_minimal() +
  theme(legend.position = "none", text = element_text(size = 21,face = "bold")) +  
  ylim(-5, 9) + 
  labs(x = "M", y = "Effects")


