# Weighted variance
weightedVariance <- function(x, w) {
  w <- w[i <- !is.na(x)]
  x <- x[i]
  sum.w <- sum(w)
  out = (sum(w*x^2) * sum.w - sum(w*x)^2) / (sum.w^2 - sum(w^2))
  out = abs(out)
  return(out)
}

# Standardized difference
Standardized_diff <- function(DM, weight, X_index, Z_index, estimand) {
  # weight: a vector of weight for each sample
  #library(SciencesPo)
  library(Matching)
  W_X = cbind(weight, DM[,X_index])
  Treated = W_X[which(DM[,Z_index]==1), ]
  Control = W_X[which(DM[,Z_index]==0), ]
  Out_before = matrix(nrow = (ncol(W_X) - 1), ncol = 5, byrow = T)
  #s = seq(from = I1, to = I2)
  for (i in 1:(ncol(W_X) - 1)) {
    Out_before[i,1] = mean(Control[,i+1])
    Out_before[i,2] = sd(Control[,i+1])
    Out_before[i,3] = mean(Treated[,i+1])
    Out_before[i,4] = sd(Treated[,i+1])
    if (sd(W_X[,i+1]) != 0) {
      Out_before[i,5] = (abs(mean(Control[,i+1]) - mean(Treated[,i+1]))/sd(W_X[,i+1]))*100
    } else if (sd(W_X[,i+1]) == 0) {
      Out_before[i,5] = NA
    }
  }
  colnames(Out_before) <- c("Control_mean", "Control_sd", "Case_mean", "Case_sd", "Before_S/D")
  if (estimand == "ATE") {
    Out_after = ATE_Weight_SD(Treated, Control, 2, ncol(W_X), 1)
  } else if (estimand == "ATT") {
    Out_after = ATT_Weight_SD(Treated, Control, 2, ncol(W_X), 1)
  } else if (estimand == "ATC") {
    Out_after0 = ATT_Weight_SD(case = Control, control = Treated, 2, ncol(W_X), 1)
    Out_after = cbind(Out_after0[,3:4], Out_after0[,1:2], Out_after0[,5])
    colnames(Out_after) = c("control.mean", "control.sd", "treated.mean", "treated.sd", "weighted.S/D")
  }
  Out = cbind(Out_before, Out_after[,5], Out_after[,c(1,3)])
  colnames(Out) = c("Control_mean", "Control_sd", "Case_mean",
                    "Case_sd", "Before_S/D", "Weighted_S/D", "Weight_control_mean","Weight_case_mean")
  return(Standerdized_diff = Out)
}

ATE_Weight_SD <- function(case, control, I1, I2, W_I) {
  #library(SciencesPo)
  s = seq(from = I1, to = I2)
  out = matrix(nrow = (I2-I1+1), ncol = 4, byrow = T)
  for (i in 1:(I2-I1+1)) {
    out[i,1] = sum(control[,s[i]]*control[,W_I])/sum(control[,W_I])
    out[i,2] = sqrt(weightedVariance(control[,s[i]], control[,W_I]))
    out[i,3] = sum(case[,s[i]]*case[,W_I])/sum(case[,W_I])
    out[i,4] = sqrt(weightedVariance(case[,s[i]], case[,W_I]))
  }
  Std_Diff = c()
  for (i in 1:nrow(out)) {
    Std_Diff[i] = (out[i,1] - out[i,3])/sqrt(((out[i,2]^2) + (out[i,4]^2))/2)
  }
  Std_Diff = Std_Diff*100
  out = cbind(out, Std_Diff)
  colnames(out) = c("control.mean", "control.sd", "treated.mean", "treated.sd", "weighted.S/D")
  return(out)
}

ATT_Weight_SD <- function(case, control, I1, I2, W_I) {
  #library(SciencesPo)
  s = seq(from = I1, to = I2)
  out = matrix(nrow = (I2-I1+1), ncol = 4, byrow = T)
  for (i in 1:(I2-I1+1)) {
    out[i,1] = sum(control[,s[i]]*control[,W_I])/sum(control[,W_I])
    out[i,2] = sqrt(weightedVariance(control[,s[i]], control[,W_I]))
    out[i,3] =  mean(case[,s[i]], na.rm = TRUE)
    out[i,4] = sd(case[,s[i]], na.rm = TRUE)
  }
  Std_Diff = c()
  for (i in 1:nrow(out)) {
    Std_Diff[i] = (out[i,1] - out[i,3])/out[i,2]
  }
  Std_Diff = Std_Diff*100
  out = cbind(out, Std_Diff)
  colnames(out) = c("control.mean", "control.sd", "treated.mean", "treated.sd", "weighted.S/D")
  return(out)
}

ATE_Weight_diff <- function(DM, weight, X_index, Z_index, std) {
  # std = c("std.diff", "std.norm")
  W_X = cbind(weight, DM[,X_index])
  case = W_X[which(DM[,Z_index]==1), ]
  control = W_X[which(DM[,Z_index]==0), ]
  I1 = 2
  I2 = ncol(W_X)
  W_I = 1
  s = seq(from = I1, to = I2)
  out = matrix(nrow = (I2-I1+1), ncol = 2, byrow = T)
  sd = out
  for (i in 1:(I2-I1+1)) {
    out[i,1] = sum(control[,s[i]]*control[,W_I])/sum(control[,W_I])
    out[i,2] = sum(case[,s[i]]*case[,W_I])/sum(case[,W_I])
    sd[i,1] = weightedVariance(control[,s[i]], control[,W_I])
    sd[i,2] = weightedVariance(case[,s[i]], case[,W_I])
  }
  #Std_Diff = out[,1] - out[,2]
  if (std == "std.diff") {
    Std_Diff = c()
    for (i in 1:nrow(out)) {
      Std_Diff[i] = (out[i,1] - out[i,2])/sqrt(((sd[i,1]) + (sd[i,2]))/2)
    }
    Std_Diff = Std_Diff*100
  } else if (std == "std.norm") {
    Std_Diff = (out[,1] - out[,2])*100
  }
  out = cbind(out, Std_Diff)
  colnames(out) = c("control.mean", "treated.mean", "weighted.diff")
  return(out)
}

# ATE estimation
ATE_infer= function(DM, weight, X_index, Z_index, normalize = T) {
  # always put y at the first column
  DM_weight = cbind(DM, weight)
  W_I = ncol(DM_weight)
  Treated = DM_weight[which(DM_weight[,Z_index]==1), ]
  control = DM_weight[which(DM_weight[,Z_index]==0), ]
  if(normalize == T) {
    ATE_case = sum(Treated[,1]*Treated[,W_I])/sum(Treated[,W_I])
    ATE_control = sum(control[,1]*control[,W_I])/sum(control[,W_I])
  } else if(normalize == F) {
    ATE_case = mean(Treated[,1]*Treated[,W_I])
    ATE_control = mean(control[,1]*control[,W_I])
  }
  ATE_after = ATE_case - ATE_control
  return(ATE_after)
}



ATE_est_averageATT.ATC = function(kbl_fit, X_index, Z_index, est_method){
  DM = kbl_fit[[2]]
  Z = DM[, Z_index]
  kbl_att = ATE_infer(DM, kbl_fit[[1]][,1], X_index, Z_index)
  kbl_atc = ATE_infer(DM, kbl_fit[[1]][,2], X_index, Z_index)
  kbl_ate = mean(Z)*kbl_att + (1-mean(Z))*kbl_atc
  out = c(kbl_ate, kbl_att, kbl_atc)
  return(out)
}

ATE_est_averageATT.ATC_withdata = function(DM, ATTweight, ATCweight,
                                           X_index, Z_index){
  #DM = kbl_fit[[2]]
  Z = DM[, Z_index]
  kbl_att = ATE_infer(DM, ATTweight, X_index, Z_index)
  kbl_atc = ATE_infer(DM, ATCweight, X_index, Z_index)
  kbl_ate = mean(Z)*kbl_att + (1-mean(Z))*kbl_atc
  out = c(kbl_ate, kbl_att, kbl_atc)
  return(out)
}

Find_bias_variance_rmse = function(DataVector, True_value, percent_bias) {
  if (percent_bias == TRUE) {
    bias = 100*(mean(DataVector, na.rm = T) - True_value)/True_value
  } else if (percent_bias == FALSE) {
    bias = mean(DataVector, na.rm = T) - True_value
  }
  variance = var(DataVector, na.rm = T)
  rmse = sqrt(mean((DataVector - True_value)^2, na.rm = T))
  stat_out = c(bias, rmse, variance)
  names(stat_out) = c("bias", "RMSE", "Var")
  return(out = stat_out)
}

Col_bias_variance_rmse = function(DM, True_value, percent_bias = TRUE) {
  if (length(True_value) == 1) {
    True_value_vec = rep(True_value, ncol(DM))
  } else {
    True_value_vec = True_value
  }
  out = matrix(nrow = 1, ncol = 3*ncol(DM))
  I = seq(1,3*ncol(DM), 3)
  for (i in 1:ncol(DM)) {
    out[,I[i]:(I[i]+2)] = Find_bias_variance_rmse(DM[,i], True_value_vec[i], percent_bias)
  }
  colnames(out) = rep(c("bias","RMSE","Var"), ncol(DM))
  return(out)
}

Get_convergence = function(Dlist) {
  p = c()
  for(i in 1:length(Dlist[[6]])) {
    p[i] = Dlist[[6]][[i]]$convergence
  }
  return(p)
}

# Calculate Hosmer-Lemeshow statistics for a list of data
Hosmer_Lemeshow = function(DM, no_quantile, Z_index, glm_ps, cbps_ps) {
  # k: number of parameters
  #require(generalhoslem)
  n = length(DM)
  HL_stat = matrix(nrow = n, ncol = 2, byrow = T)
  for (i in 1:n) {
    HL_stat[i,] = c(as.numeric(logitgof(DM[[i]][,Z_index], DM[[i]][,glm_ps], g = no_quantile)$stat),
                    as.numeric(logitgof(DM[[i]][,Z_index], DM[[i]][,cbps_ps], g = no_quantile)$stat))
  }
  colnames(HL_stat) = c("glm", "cbps")
  return(list(Hosmer_Lemeshow = HL_stat))
}

# Estimate average of outcome
Y_infer = function(DM, weight, Z_index, mu_est = "IPW") {
  DM_weight = cbind(DM[ ,1], DM[ ,Z_index], weight)
  if (mu_est == "HT") {
    mu_ipw = sum(DM_weight[,1]*DM_weight[,2]*DM_weight[,3])
  } else if (mu_est == "IPW") {
    mu_ipw = sum(DM_weight[,1]*DM_weight[,2]*DM_weight[,3])/sum(DM_weight[,2]*DM_weight[,3])
  }
  return(mu_ipw)
}
# Multiple Y
Y_infer_p = function(Y, Z, weight, mu_est = "IPW") {
  p = ncol(Y)
  ATE = c()
  for(i in 1:p) {
    DM = cbind(Y[,i], Z)
    ATE[i] = Y_infer(DM, weight = weight, Z_index = 2, mu_est = mu_est)
  }
  return(ATE)
}
# ATE estimation

ATE_infer_pY = function(DM, weight, I = seq(6,11)) {
  ATE = c()
  for(i in 1:length(I)) {
    new.data = cbind(DM[,I[i]], DM[,1:5])
    ATE[i] = ATE_infer(new.data, weight, 3:6, 2)
  }
  return(ATE)
}

# A specification test for the propensity score using its distribution conditional on participation

Normal_kernel = function(u) {(1/sqrt(2*pi)*exp(-(u^2)/2))}

TestforPS_CDFonParticipation <- function(DM, Kernel_function, c, PS_hat, Z) {
  n = nrow(DM)
  h = c*(n^(-1/8))
  Q = PS_hat
  e = DM[,Z] - Q
  V = 0
  for (i in 1:n) {
    for (j in 1:n) {
      if (i!=j) {
        Vij = Kernel_function((Q[i] - Q[j])/h)*e[i]*e[j]
        V_ij = Vij/h
        V = V + V_ij
      }
    }
  }
  V_n = V/(n*(n-1))
  return(test_stat = V_n)
}

# Logistic regression
logistic_weight = function(DM, Formula_fit, X_index, Z_index) {
  form = formula(Formula_fit)
  X = data.frame(DM[,X_index])
  Z = data.frame(DM[,Z_index])
  DM_XZ = data.frame(DM[,c(Z_index, X_index)])
  glm = glm(form, data = DM_XZ, family = binomial(link = "logit"))
  propensity_score = glm$fitted.values
  IPW = Z*(1/propensity_score) + (1-Z)*(1/(1-propensity_score))
  ATT_weight = Z+(1-Z)*(propensity_score/(1-propensity_score))
  ATC_weight = Z*((1-propensity_score)/propensity_score) + (1-Z)
  Out_weight = cbind(IPW, ATT_weight, ATC_weight, propensity_score)
  colnames(Out_weight) = c("IPW", "ATT_weight", "ATC_weight", "PS")
  Coef = as.numeric(glm$coefficients)
  #return(list(weight = Out_weight, data = DM, coefficient = Coef))
  return(list(weight = Out_weight, coefficient = Coef))
}

# CBPS
CBPS_weight = function(DM, Formula_fit, X_index, Z_index) {
  require(CBPS)
  DM_XZ = data.frame(DM[,c(Z_index, X_index)])
  cbps1_ATE = CBPS(Formula_fit, data = DM_XZ, ATT = 0, method = "exact", standardize = FALSE)
  cbps2_ATE = CBPS(Formula_fit, data = DM_XZ, ATT = 0, method = "over", standardize = FALSE)
  cbps1_ATT = CBPS(Formula_fit, data = DM_XZ, ATT = 1, method = "exact", standardize = FALSE)
  cbps2_ATT = CBPS(Formula_fit, data = DM_XZ, ATT = 1, method = "over", standardize = FALSE)
  cbps1_ATC = CBPS(Formula_fit, data = DM_XZ, ATT = 2, method = "exact", standardize = FALSE)
  cbps2_ATC = CBPS(Formula_fit, data = DM_XZ, ATT = 2, method = "over", standardize = FALSE)
  Coef = cbind(as.numeric(cbps1_ATE$coefficients), as.numeric(cbps1_ATT$coefficients),
               as.numeric(cbps1_ATC$coefficients), as.numeric(cbps2_ATE$coefficients),
               as.numeric(cbps2_ATT$coefficients), as.numeric(cbps2_ATC$coefficients))
  propensity_score = cbind(cbps1_ATE$fitted.values, cbps1_ATT$fitted.values,
                           cbps1_ATC$fitted.values, cbps2_ATE$fitted.values,
                           cbps2_ATT$fitted.values, cbps2_ATC$fitted.values)
  colnames(Coef) = colnames(propensity_score) = c("exact CBPS ATE", "exact CBPS ATT",
                                                  "exact CBPS ATC", "over CBPS ATE",
                                                  "over CBPS ATT", "over CBPS ATC")
  Out_weight = cbind(cbps1_ATE$weights, cbps1_ATT$weights, cbps1_ATC$weights,
                     cbps2_ATE$weights, cbps2_ATT$weights, cbps2_ATC$weights)
  colnames(Out_weight) = c("IPW", "ATT_weight", "ATC_weight",
                           "IPW_over", "ATT_weight_over", "ATC_weight_over")
  #return(list(weight = Out_weight, data_set = DM, propensity_score = propensity_score,
    #          coef = Coef))
  return(list(weight = Out_weight, propensity_score = propensity_score,
              coef = Coef))
}

standardize = function(DM, column = TRUE) {
  if (column == TRUE) {
    M = colMeans(DM)
    S = apply(DM, 2, sd)
    Xm = matrix(nrow = nrow(DM), ncol = ncol(DM))
    for(i in 1:ncol(DM)) {
      Xm[,i] = (DM[,i] - M[i])/S[i]
    }
  } else if (column == FALSE) {
    M = rowMeans(DM)
    S = apply(DM, 1, sd)
    Xm = matrix(nrow = nrow(DM), ncol = ncol(DM))
    for(i in 1:nrow(DM)) {
      Xm[i,] = (DM[i,] - M[i])/S[i]
    }
  }
  return(list(X = Xm, x.mean = M, x.sd = S))
}

ATE.calculate = function(Y, Z, weight) {
  ATE = c()
  for(i in 1:ncol(Y)) {
    DM = cbind(Y[,i], Z)
    ATE[i] = ATE_infer(DM, weight, 2, 2)
  }
  return(ATE)
}

CBPS_logistic = function(X, Z, Y, method) {
  if (method == "logistic") {
    logistic = glm(Z~X, family = "binomial")
    w = weight.calculate.ps(Z = Z, ps = logistic$fitted.values, standardize = TRUE)
    std = abs(Fitted.strata.diff(X = X, Z = Z, ps = logistic$fitted.values, weight = w$weight, ej))
    abs.d = abs(Fitted.strata.diff(X = X, Z = Z, ps = logistic$fitted.values, weight = w$weight, ej, std = "std.norm"))
    ATE = ATE.calculate(Y, Z, w$weight)
    ps = logistic$fitted.values
    weight = w$weight
    coef = logistic$coefficients
  } else if (method == "CBPS") {
    cbps = CBPS(Z~X, ATT = 0, method = "exact")
    #w = weight.calculate.ps(Z = Z, ps = cbps$fitted.values, standardize = TRUE)
    std = abs(Fitted.strata.diff(X = X, Z = Z, ps = cbps$fitted.values, weight = cbps$weights, ej))
    abs.d = abs(Fitted.strata.diff(X = X, Z = Z, ps = cbps$fitted.values, weight = cbps$weights, ej, std = "std.norm"))
    ATE = ATE.calculate(Y, Z, cbps$weights)
    ps = cbps$fitted.values
    weight = cbps$weights
    coef = cbps$coefficients
  }
  out = list(std.diff = std, norm.diff = abs.d, ATE = ATE, ps = ps, weight = weight, coef = coef)
  return(out)
}

SD.summary = function(SD) {
  iter = length(SD)
  p = nrow(SD[[1]])
  q = ncol(SD[[1]])
  SD.all = matrix(nrow = iter, ncol = p*q)
  for(i in 1:iter) {
    SD.all[i,] = reform.matrix1(SD[[i]])
  }
  SD.mean = colMeans(abs(SD.all), na.rm = T)
  SD.out = reform.matrix2(SD.mean, q)
  return(SD.out)
}

beta.destandardize = function(beta, mu, sigma) {
  beta = as.vector(beta)
  b0 = beta[1]
  b1 = beta[-1]
  mu = as.vector(mu)
  sigma = as.vector(sigma)
  b1.destand = b1/sigma
  b0.destand0 = b1.destand*mu
  b0.destand = b0 - sum(b0.destand0)
  beta.out = c(b0.destand, b1.destand)
  return(beta.out)
}

reform.matrix1 = function(DM) {
  p = nrow(DM)
  out = DM[1,]
  for(i in 2:p) {
    out = c(out, DM[i,])
  }
  return(out)
}

reform.matrix2 = function(DM, q) {
  p = length(DM)/q
  I1 = seq(1, length(DM), q)
  I2 = seq(q, length(DM), q)
  out = DM[I1[1]:I2[1]]
  for(i in 2:p) {
    out = rbind(out, DM[I1[i]:I2[i]])
  }
  return(out)
}

SD.summary = function(SD) {
  iter = length(SD)
  p = nrow(SD[[1]])
  q = ncol(SD[[1]])
  SD.all = matrix(nrow = iter, ncol = p*q)
  for(i in 1:iter) {
    SD.all[i,] = reform.matrix1(SD[[i]])
  }
  SD.mean = colMeans(abs(SD.all), na.rm = T)
  SD.out = reform.matrix2(SD.mean, q)
  return(SD.out)
}

Epanechnikov_kernel_parameter = function(p,q) {
  a = (-4)/((q-p)^2)
  b= 4*(q+p)/((q-p)^2)
  c = -(4*p*q)/((q-p)^2)
  #x = seq(p,q,0.001)
  #y = c + b*x + a*x^2
  #out = list(x = x, y = y, para = c(a, b, c))
  para = c(a,b,c)
  return(para)
}

Converge_plot = function(DM) {
  library(ggplot2)
  n = nrow(DM)
  newResults = as.data.frame(cbind(seq(1,n), DM))
  colnames(newResults) = c("beta","A","B","C","D","E","F1","G","H","I","J","K","L","M")
  P = ggplot(data = newResults, aes(x = beta)) +
    geom_line(aes(y = A))  +
    geom_line(aes(y = B)) +
    geom_line(aes(y = C)) +
    geom_line(aes(y = D)) +
    geom_line(aes(y = E)) +
    geom_line(aes(y = F1)) +
    geom_line(aes(y = G)) +
    geom_line(aes(y = H)) +
    geom_line(aes(y = I)) +
    geom_line(aes(y = J)) +
    geom_line(aes(y = K)) +
    geom_line(aes(y = L)) +
    geom_line(aes(y = M)) +
    labs(y = "X")
  return(P)
}

# Function to calculate the overall and strata standardized difference statistics
Standardized.Diff.Stat = function(std, ej, p1 = 15, p2 = 5) {
  std = abs(std)
  std.w0 = std[-1,]
  if (sum(is.na(std[1,])) == length(std[1,])) {
    all.max = NA
    all.mean = NA
  } else {
    all.max = max(std[1,], na.rm = T)
    all.mean = mean(std[1,], na.rm = T)
  }
  if (sum(is.na(std[-1,])) == nrow(std.w0)*ncol(std.w0)) {
    strata.mean = NA
    strata.max = NA
    strata.w.mean = NA
  } else {
    # Calculate weighted strata mean
    penalty = c(p1, rep(p2, (nrow(ej) - 2)), p1)
    std.w = std.w0 - penalty
    std.w[std.w < 0] = 0
    strata.mean = mean(std.w0, na.rm = T)
    strata.max = max(std.w0, na.rm = T)
    strata.w.mean = mean(std.w, na.rm = T)
  }
  out = c(all.max = all.max, strata.mean = strata.mean, strata.max = strata.max,
          strata.w.mean = strata.w.mean, all.mean = all.mean)
  return(out)
}

Standardized.Diff.StrataStat = function(std, ej, p1 = 15, p2 = 5) {
  std = abs(std)
  std.w0 = std[-1,]
  if (sum(is.na(std[-1,])) == nrow(std.w0)*ncol(std.w0)) {
    strata.mean = NA
    strata.max = NA
    strata.median = NA
    strata.w.mean = NA
  } else if (sum(is.na(std[1,])) == length(std[1,])) {
    overall.mean = NA
  } else {
    # Calculate weighted strata mean
    penalty = c(p1, rep(p2, (nrow(ej) - 2)), p1)
    std.w = std.w0 - penalty
    std.w[std.w < 0] = 0
    strata.mean = mean(std.w0, na.rm = T)
    strata.max = max(std.w0, na.rm = T)
    strata.median = median(std.w0, na.rm = T)
    strata.w.mean = mean(std.w, na.rm = T)
    overall.mean = mean(std[1,], na.rm = T)
  }
  out = c(strata.mean = strata.mean, strata.max = strata.max, strata.median = strata.median,
          strata.w.mean = strata.w.mean, overall.mean = overall.mean)
  return(out)
}

Weight_stability_stat = function(w) {
  if(sum(is.na(w)) == length(w)) {
    out = rep(NA, 3)
  } else {
    mean_w = mean(w, na.rm = TRUE)
    sd_w = sd(w, na.rm = TRUE)
    quantile_w = quantile(w, probs = c(0.95, 0.99), na.rm = TRUE)
    out = c((sd_w/mean_w), quantile_w)
  }
  names(out) = c("cv","p95","p99")
  return(out)
}

# Function to generate interaction terms in a given dataset
Generate_interaction = function(DM, x, m = 2) {
  combine_number = combn(x,m)
  n = nrow(DM)
  p = ncol(combine_number)
  DM_inter = matrix(nrow = n, ncol = p, byrow = F)
  for(i in 1:p) {
    DM_inter[,i] = DM[,combine_number[1,i]]*DM[,combine_number[2,i]]
  }
  return(DM_interaction = DM_inter)
}

# Function to calculate Komogrov statistics
Komogrov_calculate = function(X, Z, weight) {
  X0 = X
  X = cbind(rep(1,nrow(X)), X)
  control = X0[which(Z==0),]
  case = X0[which(Z==1),]
  w_control = weight[which(Z==0)]
  w_case = weight[which(Z==1)]
  Komogrov = c()
  for(i in 1:ncol(X0)) {
    Komogrov[i] = as.numeric(ks.test(control[,i]*w_control, case[,i]*w_case, exact = FALSE)$statistic)
  }
  return(Komogrov)
}

# Function to calculate standardized difference in each given strata of data
Fitted_diff = function(X, Z, beta, ej) {
  DM = cbind(Z, X)
  #X = cbind(rep(1, nrow(X)), X)
  #fit_t = exp(X%*%matrix(beta_t, ncol = 1))/(1+exp(X%*%matrix(beta_t, ncol = 1)))
  fit_value = exp(X%*%matrix(beta, ncol = 1))/(1+exp(X%*%matrix(beta, ncol = 1)))
  weight = Z/fit_value + (1-Z)/(1-fit_value)
  SD = matrix(nrow = nrow(ej), ncol = ncol(X)-1, byrow = T)
  for(i in 1:nrow(ej)) {
    I = which(fit_value>=ej[i,1]&fit_value<ej[i,2])
    if(length(I) <= 5) {
      SD[i,] = rep(NA, ncol = ncol(X)-1)
    } else if (sum(DM[I,1]==0) <=2 | sum(DM[I,1]==1) <= 2) {
      SD[i,] = rep(NA, ncol = ncol(X)-1)
    } else if (sum(DM[I,1]==0) == length(DM[I,1]) | sum(DM[I,1]==1) == length(DM[I,1])) {
      SD[i,] = rep(NA, ncol = ncol(X)-1)
    } else {
      SD[i,] = Standardized_diff(DM[I,], weight[I],2:ncol(DM),1,"ATE")[,6]
    }
  }
  SD_total =  Standardized_diff(DM, weight,2:ncol(DM),1,"ATE")[,6]
  out = rbind(SD_total, SD)
  return(out)
}
#' Calculate Local Covariate Balance in Pre-specified Local Neighborhoods
#'
#' This function evaluates the global and local balance by standardized difference
#' or non-standardized difference of covariates after inverse probability weighting (IPW) weighting. The IPW is used for ATE estimation.
#' The global balance is the difference measures calculated in all the subjects.
#' The local balance is the difference measures calculated in each pre-defined local neighborhood.
#'
#' @param X The covariate matrix with the rows corresponding to the subjects/participants and the columns corresponding to the covariates.
#' @param Z The binary treatment indicator. A vector with 2 unique numeric values in 0 = untreated and 1 = treated.
#'          Note each subject/participant must has a value of the treatment indicator,
#'          thus its length should be the same as the number of all subjects/participants.
#' @param ps The propensity score used to determine the propensity score (PS) stratified subpopulation in each local neighborhoods.
#'           A vector with numeric values greater than 0 and less than 1 with the legnth of the number of all subjects/participants.
#' @param weight The IPW weight used for estimating average treatement effect (ATE).
#'               It is calculated as Z/ps + (1-Z)/(1-ps), where Z is the treatment indicator and ps is the propensity score.
#'               A vector with positive values and the legnth of the number of all subjects/participants.
#'               If weight = NULL, the IPW for ATE will be calculated using the input propensity score. The weight will be normalized
#'               to sum to 1 within treated/untreated group. weight = NULL by default.
#' @param ej The matrix of the local neighborhoods, which contains two columns of positive values greater or equal to 0 and less or equal to 1.
#'           The rows of ej represent the neighborhoods. The first column is the start point of the local neighborhoods.
#'           The second column is the end point of the neighborhoods.
#' @param std The difference measure that evaluating the covariate balance. The default measure "std.diff" is the
#'            standardized difference of covariates, which is the difference in weighted means of the two treatment groups,
#'            divided by the pooled weighted standard deviation expressed as a percentage. The measure "std.norm" is the
#'            difference in weighted means of the two treatment groups. "std" can only take values of "std.diff" or "std.norm".
#'
#' @return The matrix whose each column representing the difference measure for each covariates.
#'         The first row represents the global balance. The rest of the row represents the local
#'         balance corresponding to the local neighborhoods. The global or local balance
#'         are standardized difference of each covariate in the whole or propensity score stratified (PS-stratified)
#'         sub-populations, and expressed as a percentage.
#'
#' @references Li, L. and Greene, T. (2013) A weighting analogue to pair matching in propensity score analysis.
#'             The International Journal of Biostatistics, 9, 215-234.
#'
#' @examples
#' KS = Kang_Schafer_Simulation(n = 2000, seeds = 5050)
#' # Define local neighborhoods
#' ej = cbind(seq(0,0.7,0.1),seq(0.3,1,0.1))
#' print(ej)
#' # Misspecified propensity score model
#' require(CBPS)
#' Z = KS$Data[,2]
#' X = KS$Data[,7:10]
#' cbps.fit = CBPS(Z~X, ATT = 0, method = "exact")
#' # Global and local balance of CBPS estimated propensity score
#' Fitted.strata.diff(X = X, Z = Z, ps = cbps.fit$fitted.values,
#'                    weight = cbps.fit$weights, ej = ej, std = "std.diff")
#' # Global and local balance of the true propensity score
#' true.ps = KS$Data[,11]
#' Fitted.strata.diff(X = X, Z = Z, ps = true.ps, ej = ej, std = "std.diff")
#'
#' @export
Fitted.strata.diff = function(X, Z, ps, weight = NULL, ej, std = "std.diff") {
  if (std != "std.diff" & std != "std.norm") {
    stop("std can only take values of std.diff or std.norm")
  }
  DM = cbind(Z, X)
  #X = cbind(rep(1, nrow(X)), X)
  fit_value = ps
  if (is.null(weight)) {
    weight = weight.calculate.ps(Z = Z, ps = ps, standardize = TRUE)$weight
  } else if (!is.null(weight)) {
    weight = weight
  }
  #weight = Z/fit_value + (1-Z)/(1-fit_value)
  SD = matrix(nrow = nrow(ej), ncol = ncol(X), byrow = T)
  for(i in 1:nrow(ej)) {
    I = which(fit_value>=ej[i,1]&fit_value<ej[i,2])
    if(length(I) <= 5) {
      SD[i,] = rep(NA, ncol = ncol(X)-1)
    } else if (sum(DM[I,1]==0)<=2 | sum(DM[I,1]==1)<=2) {
      SD[i,] = rep(NA, ncol = ncol(X)-1)
    } else if (sum(DM[I,1]==0) == length(DM[I,1]) | sum(DM[I,1]==1) == length(DM[I,1])) {
      SD[i,] = rep(NA, ncol = ncol(X)-1)
    } else {
      if(std == "std.diff") {
        SD[i,] = Standardized_diff(DM[I,], weight[I],2:ncol(DM),1,"ATE")[,6]
      } else if(std == "std.norm") {
        SD[i,] = ATE_Weight_diff(DM[I,], weight[I], 2:ncol(DM), 1, std = "std.norm")[,3]
      }
    }
  }
  if(std == "std.diff") {
    SD_total =  Standardized_diff(DM, weight,2:ncol(DM),1,"ATE")[,6]
  } else if(std == "std.norm") {
    SD_total = ATE_Weight_diff(DM, weight,2:ncol(DM),1,std = "std.norm")[,3]
  }
  out = rbind(SD_total, SD)
  X.name = colnames(X)
  if (is.null(X.name)) {
    X.name = paste0("X",seq(1,ncol(X)))
  } else {
    X.name = X.name
  }
  S.name = rownames(ej)
  if (is.null(S.name)) {
    S.name = c()
    for (i in 1:nrow(ej)) {
      S.name[i] = paste0(ej[i,1], "-", ej[i,2])
    }
  } else {
    S.name = S.name
  }
  S.name = c("global_bal", S.name)
  rownames(out) = S.name
  colnames(out) = X.name
  return(out)
}

# Function to calculate standardized difference in each given strata of data
weight.calculate = function(X, Z, beta, standardize = FALSE, Intercept = TRUE) {
  beta = matrix(beta, ncol = 1)
  n = nrow(X)
  if (Intercept == TRUE) {
    X = cbind(rep(1,n), X)
  } else if (Intercept == FALSE) {
    X = X
  }
  ps = exp(X%*%beta)/(1+exp(X%*%beta))
  ps = as.vector(ps)
  w = Z/ps + (1-Z)/(1-ps)
  if (standardize == TRUE) {
    I0 = which(Z == 0)
    I1 = which(Z == 1)
    weight = w
    weight[I0] = w[I0]/sum(w[I0])
    weight[I1] = w[I1]/sum(w[I1])
  } else if (standardize == FALSE) {
    weight = w
  }
  return(list(weight = weight, ps = ps))
}

#' The Inverse Probability Weight (IPW) for Average Treatment Effect (ATE)
#'
#' This function calculates the IPW for ATE as Tr/ps + (1-Tr)/(1-ps), where Tr is the treatment
#' indicator and ps is the propensity score.
#'
#' @param Z The binary treatment indicator. A vector with 2 unique numeric values in 0 = untreated and 1 = treated.
#' @param ps The propensity score used to calculate the IPW weight. A vector with numeric values greater than 0 and less than 1.
#'           The length of "ps" should equal to the length of "Z".
#' @param standardize If standardize the calculated IPW or not. When standardize = FALSE, the weights are calculated
#'                    by the IPW definition shown in the description. standardize = TRUE by default, which makes the
#'                    the sum of weights be 1 in the untreated/treated group.
#'
#' @return A list containing the following components:
#'         \itemize{
#'           \item{"weight"}: The IPW weights calculated by the input propensity score
#'           \item{"ps"}: The input propensity score
#'         }
#'
#' @examples
#' KS = Kang_Schafer_Simulation(n = 1000, seeds = 5050)
#' tr = KS$Data[,2]
#' true.ps = KS$Data[,11]
#' true.w = weight.calculate.ps(Z = tr, ps = true.ps, standardize = TRUE)
#' summary(true.w$weight)
#' c(sum(true.w$weight[which(tr==0)]), sum(true.w$weight[which(tr==1)]))
#'
#' @export

weight.calculate.ps = function(Z, ps, standardize = TRUE) {
  n = length(ps)
  probs.min = 1e-6
  ps = pmin(1-probs.min, ps)
  ps = pmax(probs.min, ps)
  ps = as.vector(ps)
  w = Z/ps + (1-Z)/(1-ps)
  if (standardize == TRUE) {
    I0 = which(Z == 0)
    I1 = which(Z == 1)
    weight = w
    weight[I0] = w[I0]/sum(w[I0])
    weight[I1] = w[I1]/sum(w[I1])
  } else if (standardize == FALSE) {
    weight = w
  }
  return(list(weight = weight, ps = ps))
}

# SD
SD.organize = function(SD) {
  iter = length(SD)
  k = dim(SD[[1]])[1] - 1
  p = dim(SD[[1]])[2]
  cov.name = paste0("cov",seq(1,p))
  strata.name = c("all",paste0("K",seq(1,k)))
  out = list()
  for(i in 1:p) {
    out[[i]] = matrix(nrow = iter, ncol = (k + 1))
    for(j in 1:iter) {
      out[[i]][j,] = SD[[j]][,i]
    }
    colnames(out[[i]]) = strata.name
  }
  names(out) = cov.name
  return(out)
}



