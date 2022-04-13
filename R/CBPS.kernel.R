#' Propensity Score Analysis with Local Balance 1 (PSLB 1) Estimation
#'
#' PSLB1 estimates propensity scores by implementing a flexible form of the covariate balancing propensity score (CBPS) using kernel PCA,
#' and tunes parameters (the bandwidth of the Gaussian kernel and the number of PCs) by the tuning for local balance (T4LB) algorithm, which finds
#' the model with the best local balance (minimum absolute standardized difference (S/D) in the PS-stratified sub-populations) among the model pool. The method
#' searches for the propensity score model with the best local balance while controling the global balance. The estimation is considered as "fail" if the minimum
#' absolute S/D of input covariates in the whole population is more than 10% (uncontrolled global balance), and the function output NA
#' value with an error message. The local balance is evaluated by a statistic (mean or max) of the absolute S/D in the PS-stratified sub-populations.
#' The method only takes binary treatment.
#'
#' @param X The covariate matrix with the rows corresponding to the subjects/participants and the columns corresponding to the covariates.
#'          The covariate matrix X can not contain missing value. The function will stop if NA value is detected in X.
#' @param Z The binary treatment indicator. A vector with 2 unique numeric values in 0 = untreated and 1 = treated.
#'          Its length should be the same as the number of all subjects/participants without missing value.
#'          The function will stop if NA value is detected in Z.
#' @param n_sigma The number of bandwidth value of the Gaussian kernel function used to generate the feature space.
#'                See details in the description of function \code{\link{Gaussian_Kernel_feature}}. Default to 20.
#' @param ej The matrix of the local neighborhoods with its rows representing the neighborhoods. It contains two columns of values greater or equal to 0 and less or equal to 1.
#'           The first columns are the start and the second column are the end point of the local neighborhoods.
#' @param selectInX The matrix used to evaluate the local balance. selectInX = "cov" by default, which evaluates the local covariate balance on the
#'                  input covariate matrix X. If selectInX = "feature", the local covariate balance is evaluated on the feature space transformed by kernel PCA.
#' @param k The number of top K eigen value calculated. Note that k should be smaller than the sample size. The minimum of k and the sample size will be used.
#'          Defaults to 500. See details in the description of function \code{\link{Gaussian_Kernel_feature}}.
#' @param method Choose "exact" to fit the justidentified CBPS model; choose "over" to fit the overidentified CBPS model. See details in Imai and Ratkovic (2014).
#'               Default is "exact".
#' @param standardize Default is TRUE, which normalizes weights to sum to 1 within treated/untreated group.
#'                    Set to FALSE to return IPW weights for ATE.
#' @param criteria Choose "mean" to use the mean of the absolute S/D among all covariates and neighborhoods for the tuning parameters by T4LB;
#'                 choose "max" to use the max of the absolute S/D among all covariates and neighborhoods for the tuning parameters by T4LB.
#'                 Default is "mean".
#'
#' @return A list containing the following components:
#'         \itemize{
#'           \item{"weight"}: The IPW weights for binary ATE calculated by the fitted propensity score given by Tr/ps + (1-Tr)/(1-ps),
#'                            where Tr is the treatment assignment and ps is the fitted propensity score. This expression for weight
#'                            is before standardization (i.e. with standardize = FALSE). Standardization will make weights sum to 1
#'                            within untreated/treated group.
#'           \item{"propensity.score"}: The fitted propensity score
#'           \item{"balance"}: The matrix containing the global and local balance of the estimated propensity score. The first row contains the absolute S/D of each covariates
#'                             in the whole study population, which represents the global covariate balance. The rest rows contain the absolute
#'                             S/D of each covariates in the PS-stratified sub-population corresponding to the loal neighborhoods of "ej", which
#'                             represent the local covariate balance.
#'           \item{"coefficients"}: The coefficient of the fitted propensity score model
#'           \item{"feature"}: The selected features by kernel PCA with the best loal covariate balance.
#'           \item{"best.para"}: The parameters (bandwidth for Gaussian kernel and number of PCs) chosen with the best local covariate balance.
#'           \item{"treat"}: The treatment assignment vector used
#'         }
#'
#' @references Imai, K. and Ratkovic, M. (2014) Covariate balancing propensity score. Journal of the
#'             Royal Statistical Society: Series B (Statistical Methodology), 76, 243-263.
#'
#' @examples
#' KS = Kang_Schafer_Simulation(n = 1500, seeds = 5050)
#' # Misspecified propensity score model
#' X = KS$Data[,7:10]
#' Z = KS$Data[,2]
#' # Specify the local neighborhoods
#' ej = cbind(seq(0,0.8,0.2),seq(0.2,1,0.2))
#' print(ej)
#' # PSLB 1 fitting
#' fit = PSLB1(X = X, Z = Z, n_sigma = 10, ej = ej)
#' print(fit$balance)
#'
#' @export
PSLB1 = function(X, Z, n_sigma = 20, ej, selectInX = "cov", k = 500, method = "exact", standardize = TRUE,
                 criteria = "mean") {
  if (criteria != "mean" && criteria != "max") {
    stop("Please choose an available local balance statistic from mean and max.")
  }
  na.x = sum(is.na(X))
  na.z = sum(is.na(Z))
  if (na.x != 0) {
    stop("The covariate matrix can not contain missing value.")
  }
  if (na.z != 0) {
    stop("The treatment indicator can not contain missing value.")
  }
  fit = NA
  tryCatch({fit = CBPS.kernel(X = X, Z = Z, n_sigma = n_sigma, ej = ej, selectInX = selectInX, k = k,
                              standardize = standardize, method = method)},
           error = function(e) {
             print("PSLB1 failed due to fitting error")
           })
  if (is.na(fit[1])) {
    out = NULL
  } else if (sum(is.na(fit$exact$best.fit$best.coef$mean)) != 0) {
    print("PSLB1 failed due to uncontrolledly large global covariate balance")
    out = NULL
  } else if (sum(is.na(fit$exact$best.fit$best.coef$mean)) == 0) {
    if (criteria == "mean") {
      weight = fit$exact$best.fit$best.weight[,1]
      ps = fit$exact$best.fit$best.ps[,1]
      coef = fit$exact$best.fit$best.coef$mean
      best.para0 = fit$exact$best.index$parameter[fit$exact$best.index$best.index[1],]
    } else if (criteria == "max") {
      weight = fit$exact$best.fit$best.weight[,2]
      ps = fit$exact$best.fit$best.ps[,2]
      coef = fit$exact$best.fit$best.coef$max
      best.para0 = fit$exact$best.index$parameter[fit$exact$best.index$best.index[2],]
    }
    # best parameters selected
    best.para = cbind(best.para0[1:2],c(fit$sigma[best.para0[1]], best.para0[3]))
    colnames(best.para) = c("index","value")
    rownames(best.para) = c("sigma", "l")
    # feature space
    I = fit$r[best.para0[1], best.para0[2]]
    feature = fit$feature[[best.para0[1]]][,1:I]
    # Global and local balance
    std = abs(Fitted.strata.diff(X = X, Z = Z, ps = ps, weight = weight, ej = ej))
    out = list(weight = weight, propensity.score = ps, balance = std, coefficients = coef,
               feature = feature, best.para = best.para, treat = Z)
  }
  return(out)
}
# PSLB 1 fitting with or without cross validation
# this function can include other covariate transformations as cov_input
CBPS.kernel = function(X, Z, n_sigma, n_fold, ej, selectInX, k = 500, cov_input = NULL, CrossVal = FALSE,
                       p1 = 15, p2 = 5, standardize = TRUE, method = "exact") {
  set.seed(5050)
  require(CBPS)
  kernel.fit = Gaussian_Kernel_feature(X, n_sigma = n_sigma, k = k)
  if (selectInX == "cov") {
    X.std = X
  } else if (selectInX == "feature") {
    X.std = NULL
  } else if (!is.null(cov_input) && selectInX == "cov_input") {
    X.std = cbind(X, cov_input)
  }
  if (is.null(cov_input) == TRUE) {
    feature.list = kernel.fit$features
    r = kernel.fit$r
    number_input = 0
  } else if (is.null(cov_input) == FALSE) {
    cov_input_sd = standardize(cov_input)$X
    feature.list = list()
    for(i in 1:length(kernel.fit$features)) {
      feature.list[[i]] = list()
      for(j in 1:length(kernel.fit$features[[i]])) {
        feature.list[[i]][[j]] = cbind(kernel.fit$features[[i]][[j]], cov_input_sd)
      }
    }
    number_input = dim(cov_input_sd)[2]
    r = kernel.fit$r + number_input
  }
  if (method == "exact") {
    exact.fit = CBPS.kernel.CV.ChoosePara(feature.list, Z, r, CrossVal = CrossVal, n_fold,
                                          method = "exact", ej, p1, p2, X.std = X.std, standardize = standardize)
  } else if (method == "over") {
    exact.fit = CBPS.kernel.CV.ChoosePara(feature.list, Z, r, CrossVal = CrossVal, n_fold,
                                          method = "over", ej, p1, p2, X.std = X.std, standardize = standardize)
  }
  feature.out = list()
  for(i in 1:length(feature.list)) {
    feature.out[[i]] = feature.list[[i]][[3]]
  }
  sigma = kernel.fit$sigma_seq
  out = list(exact = exact.fit, r = r, feature = feature.out,
             number_input = number_input, method = method, sigma = sigma)
  return(out)
}

# Function to implement Local-Cov-Balance algorithm in kernelized CBPS
CBPS.kernel.CV.ChoosePara = function(feature.list, Z, r, CrossVal, n_fold,
                                     method, ej, p1, p2, X.std, standardize = TRUE) {
  n = nrow(feature.list[[1]][[1]])
  p = rep(1,ncol(r))
  r.all = r[1,]
  for(i in 2:nrow(r)) {
    p = c(p, rep(i,ncol(r)))
    r.all = c(r.all, r[i,])
  }
  #q = rep(seq(1,5),nrow(r))
  q = rep(seq(1,3),nrow(r))
  para = cbind(p, q, r.all)
  colnames(para) = c("sigma","r","r.no")
  if (CrossVal == TRUE) {
    std.list = list()
    for(i in 1:nrow(para)) {
      if (is.null(X.std)) {
        cbps.fit0 = CBPS.kernel.fit(feature.list[[para[i,1]]][[para[i,2]]], Z, method,
                                    ej, p1, p2, X.std = feature.list[[para[i,1]]][[para[i,2]]],
                                    standardize = standardize)
      } else {
        cbps.fit0 = CBPS.kernel.fit(feature.list[[para[i,1]]][[para[i,2]]], Z, method,
                                    ej, p1, p2, X.std = X.std, standardize = standardize)
      }
      std.list[[i]] = cbps.fit0$std.summary
    }
    best.I = CBPS.BestIndex(std.list = std.list, para = para)
    best.fit = CBPS.BestFitCV(feature.list, Z, I = best.I$best.index,
                              para = para, method = method, ej = ej,
                              p1 = p1, p2 = p2, standardize = standardize)
  } else if (CrossVal == FALSE) {
    std.list = list()
    fit.list = list()
    for(i in 1:nrow(para)) {
      if (is.null(X.std)) {
        cbps.fit0 = CBPS.kernel.fit(feature.list[[para[i,1]]][[para[i,2]]], Z, method,
                                    ej, p1, p2, X.std = X.std, standardize = standardize)
      } else {
        cbps.fit0 = CBPS.kernel.fit(feature.list[[para[i,1]]][[para[i,2]]], Z, method,
                                    ej, p1, p2, X.std = X.std, standardize = standardize)
      }
      std.list[[i]] = cbps.fit0$std.summary
      fit.list[[i]] = cbps.fit0
    }
    best.I = CBPS.BestIndex(std.list = std.list, para = para)
    best.fit = CBPS.BestFitnoCV(fit.list, I = best.I$best.index, n)
  }
  out = list(best.fit = best.fit, best.index = best.I)
  return(out)
}

# kernelized CBPS
CBPS.kernel.fit = function(features, Z, method, ej, p1, p2, X.std,
                           X = NULL, standardize = TRUE) {
  # method = c("exact", "over")
  DM = data.frame(cbind(Z, features))
  r = ncol(features)
  N = paste0("X", seq(1,r))
  colnames(DM) = c("Z", N)
  form = Form.formula(N)
  form = formula(form)
  cbps_ATE = NA
  tryCatch({cbps_ATE = CBPS(form, data = DM, ATT = 0, method = method)},
            error = function(e) {
              print("CBPS fit error")
            })
  if (is.na(cbps_ATE[1])) {
    cbps.coef = rep(NA, r)
    cbps.ps = cbind(rep(NA, length(Z)), rep(NA, length(Z)))
    data = cbind(Z, features)
    colnames(cbps.ps) = c("ps","weight")
    if (is.null(X.std)) {
      std = matrix(rep(NA,(nrow(ej)+1)*ncol(feature)), nrow = (nrow(ej)+1), ncol = ncol(feature))
    } else {
      std = matrix(rep(NA,(nrow(ej)+1)*ncol(X.std)), nrow = (nrow(ej)+1), ncol = ncol(X.std))
    }
    std.summary = rep(NA, 5)
    names(std.summary) = c("all.max","strata.mean","strata.max","strata.w.mean","all.mean")
  } else {
    cbps.coef = cbps_ATE$coefficients
    if (is.null(X) == TRUE) {
      cbps.ps = cbind(cbps_ATE$fitted.values, cbps_ATE$weights)
      data = cbind(Z, features)
    } else if (is.null(X) == FALSE) {
      w = weight.calculate(X = X[,-1], Z = X[,1],
                           cbps.coef, standardize = standardize)
      cbps.ps = cbind(w$ps, w$weight)
      data = X
    }
    colnames(cbps.ps) = c("ps","weight")
    if (is.null(X.std)) {
      std = abs(Fitted.strata.diff(data[,-1], data[,1],
                                   cbps.ps[,1], cbps.ps[,2], ej))
    } else {
      std = abs(Fitted.strata.diff(X.std, data[,1],
                                   cbps.ps[,1], cbps.ps[,2], ej))
    }
    std.summary = Standardized.Diff.Stat(std = std, ej = ej, p1 = p1, p2 = p2)
  }

  cbps.fit = list(coefficient = cbps.coef, propensity.score = cbps.ps[,1],
                  weight = cbps.ps[,2], std.summary = std.summary,
                  std = std, data = data)
  return(cbps.fit)
}

# Function implement cross validation in CBPS-Kernel
CBPS.kernel.fit.CV = function(features, Z, n_fold, method,
                              ej, p1, p2, standardize = TRUE) {
  n = nrow(features)
  I = sample(seq(1,n),n)
  p = round(n/n_fold)
  I1 = seq(1,n,p)[1:n_fold]
  I2 = I1 + p -1
  I2 = c(I2[1:(n_fold-1)], n)
  I_fit = list()
  I_test = list()
  for(i in 1:n_fold) {
    I_fit[[i]] = I[-(I1[i]:I2[i])]
    I_test[[i]] = I[I1[i]:I2[i]]
  }
  CV_fit = list()
  std.sum = matrix(nrow = n_fold, ncol = 5)
  for(i in 1:n_fold) {
    X_fit = cbind(Z, features)[I_fit[[i]],]
    X_test = cbind(Z, features)[I_test[[i]],]
    CV_fit0 = CBPS.kernel.fit(features = X_fit[,-1],
                              Z = X_fit[,1], method = method,
                              X = X_test, ej = ej,
                              p1 = p1, p2 = p2, standardize = standardize)
    CV_fit[[i]] = list(coefficient = CV_fit0$coefficient, data = CV_fit0$data,
                       std = CV_fit0$std.summary)
    std.sum[i,] = CV_fit0$std.summary
  }
  std.sum.mean = colMeans(std.sum, na.rm = TRUE)
  out = list(CV_fit = CV_fit, std.summary = std.sum.mean)
  return(out)
}

CBPS.BestIndex = function(std.list, para) {
  std.summary = matrix(nrow = nrow(para), ncol = 5)
  for(i in 1:nrow(para)) {
    std.summary[i,] = std.list[[i]]
  }
  colnames(std.summary) = c("all.max","strata.mean","strata.max","strata.w.mean","all.mean")
  if (length(which(std.summary[,1]<=10)) == 0) {
    mean.m = I1 = NA
    max.m = I2 = NA
    w.mean.m = min(std.summary[,4], na.rm = TRUE)
  } else if (length(which(std.summary[,1]<=10)) != 0) {
    mean.m = min(std.summary[which(std.summary[,1]<=10),2], na.rm = TRUE)
    max.m = min(std.summary[which(std.summary[,1]<=10),3], na.rm = TRUE)
    I1 = which(std.summary[,2] == mean.m)
    I2 = which(std.summary[,3] == max.m)
    w.mean.m = min(std.summary[which(std.summary[,1]<=10),4], na.rm = TRUE)
  }
  meanOnly.m = min(std.summary[,2], na.rm = TRUE)
  overall.mean.m = min(std.summary[,5], na.rm = TRUE)
  I = c(I1, I2,
        which(std.summary[,4] == w.mean.m),
        which(std.summary[,2] == meanOnly.m),
        which(std.summary[,5] == overall.mean.m))
  out = list(best.index = I, parameter = para)
  return(out)
}

CBPS.BestFitCV = function(feature.list, Z, I, para, method, ej, p1, p2, standardize) {
  n = nrow(feature.list[[1]][[1]])
  if (is.na(I[1]) == TRUE) {
    mean.best = max.best = list(coefficient = NA, propensity.score = rep(NA, n),
                                weight = rep(NA, n), std = NA)
  } else if (is.na(I[1]) == FALSE) {
    mean.best = CBPS.kernel.fit(feature.list[[para[I[1],1]]][[para[I[1],2]]], Z,
                                method, ej, p1, p2, standardize = standardize)
    max.best = CBPS.kernel.fit(feature.list[[para[I[2],1]]][[para[I[2],2]]], Z,
                               method, ej, p1, p2, standardize = standardize)
  }
  w.mean.best = CBPS.kernel.fit(feature.list[[para[I[3],1]]][[para[I[3],2]]], Z,
                                method, ej, p1, p2, standardize = standardize)
  MeanOnly.best = CBPS.kernel.fit(feature.list[[para[I[4],1]]][[para[I[4],2]]], Z,
                                  method, ej, p1, p2, standardize = standardize)
  overall.mean.best = CBPS.kernel.fit(feature.list[[para[I[5],1]]][[para[I[5],2]]], Z,
                                      method, ej, p1, p2, standardize = standardize)
  coef = list(mean = mean.best$coefficient, max = max.best$coefficient,
              w.mean = w.mean.best$coefficient, MeanOnly = MeanOnly.best$coefficient,
              overall.mean = overall.mean.best$coefficient)
  ps = cbind(mean.best$propensity.score, max.best$propensity.score,
             w.mean.best$propensity.score, MeanOnly.best$propensity.score,
             overall.mean.best$propensity.score)
  w = cbind(mean.best$weight, max.best$weight, w.mean.best$weight,
            MeanOnly.best$weight, overall.mean.best$weight)
  std = list(mean = mean.best$std, max = max.best$std, w.mean = w.mean.best$std,
             MeanOnly = MeanOnly.best$std, overall.mean = overall.mean.best$std)
  colnames(ps) = colnames(w) = c("mean","max","w.mean","MeanOnly","overall.mean")
  out = list(best.coef = coef, best.ps = ps, best.weight = w, best.std = std)
  return(out)
}

CBPS.BestFitnoCV = function(fit.list, I, n) {
  #n = nrow(feature.list[[1]][[1]])
  if (is.na(I[1]) == TRUE) {
    mean.best = max.best = list(coefficient = NA, propensity.score = rep(NA, n),
                                weight = rep(NA, n), std = NA)
  } else if (is.na(I[1]) == FALSE) {
    mean.best = fit.list[[I[1]]]
    max.best = fit.list[[I[2]]]
  }
  w.mean.best = fit.list[[I[3]]]
  MeanOnly.best = fit.list[[I[4]]]
  overall.mean.best = fit.list[[I[5]]]
  coef = list(mean = mean.best$coefficient, max = max.best$coefficient,
              w.mean = w.mean.best$coefficient, MeanOnly = MeanOnly.best$coefficient,
              overall.mean = overall.mean.best$coefficient)
  ps = cbind(mean.best$propensity.score, max.best$propensity.score,
             w.mean.best$propensity.score, MeanOnly.best$propensity.score,
             overall.mean.best$propensity.score)
  w = cbind(mean.best$weight, max.best$weight, w.mean.best$weight,
            MeanOnly.best$weight, overall.mean.best$weight)
  std = list(mean = mean.best$std, max = max.best$std, w.mean = w.mean.best$std,
             MeanOnly = MeanOnly.best$std, overall.mean = overall.mean.best$std)
  colnames(ps) = colnames(w) = c("mean","max","w.mean","MeanOnly","overall.mean")
  out = list(best.coef = coef, best.ps = ps, best.weight = w, best.std = std)
  return(out)
}

CBPS.kernelall = function(X, Z, n_sigma, n_fold, ej,
                          p1 = 15, p2 = 5, standardize = TRUE) {
  set.seed(5050)
  require(CBPS)
  kernel.fit = Gaussian_Kernel_feature(X, n_sigma = n_sigma)
  feature.list = kernel.fit$features
  r = kernel.fit$r
  exact.CV = CBPS.kernel.CV.ChoosePara(feature.list, Z, r, CrossVal = TRUE, n_fold,
                                       method = "exact", ej, p1, p2, standardize = standardize)
  over.CV = CBPS.kernel.CV.ChoosePara(feature.list, Z, r, CrossVal = TRUE, n_fold,
                                      method = "over", ej, p1, p2, standardize = standardize)
  exact.noCV = CBPS.kernel.CV.ChoosePara(feature.list, Z, r, CrossVal = FALSE, n_fold,
                                         method = "exact", ej, p1, p2, standardize = standardize)
  over.noCV = CBPS.kernel.CV.ChoosePara(feature.list, Z, r, CrossVal = FALSE, n_fold,
                                        method = "over", ej, p1, p2, standardize = standardize)
  feature.out = list()
  for(i in 1:length(feature.list)) {
    feature.out[[i]] = feature.list[[i]][[5]]
  }
  out = list(exact.CV = exact.CV, over.CV = over.CV,
             exact.noCV = exact.noCV, over.noCV = over.noCV,
             r = r, feature = feature.out)
  return(out)
}

CBPS.kernel.best = function(cbps.k, X, Z, ej, X1) {
  p = length(cbps.k$cbps.fit)
  r = seq(1,length(cbps.k$cbps.fit[[1]]))
  sigma = rep(1,length(cbps.k$cbps.fit[[1]]))
  for(i in 2:p) {
    r = c(r, seq(1,length(cbps.k$cbps.fit[[i]])))
    sigma = c(sigma, rep(i,length(cbps.k$cbps.fit[[i]])))
  }
  std.out = list()
  all.max = c()
  strata.mean = c()
  strata.max = c()
  strata.w.mean = c()
  for(i in 1:length(r)) {
    std = Fitted.strata.diff(X, Z, cbps.k$cbps.fit[[sigma[i]]][[r[i]]]$propensity.score[,1],
                             cbps.k$cbps.fit[[sigma[i]]][[r[i]]]$propensity.score[,2], ej)
    std = abs(std)
    std.out[[i]] = std
    std.w0 = std[-1,]
    # Calculate weighted strata mean
    p1 = 15
    p2 = 5
    penalty = c(p1, rep(p2, (nrow(ej) - 2)), p1)
    std.w = std.w0 - penalty
    std.w[std.w < 0] = 0
    all.max[i] = max(std[1,], na.rm = T)
    strata.mean[i] = mean(std[-1,], na.rm = T)
    strata.max[i] = max(std[-1,], na.rm = T)
    strata.w.mean[i] = mean(std.w, na.rm = T)
  }
  I = which(all.max <= 10)
  S = c(which(strata.mean == min(strata.mean[I])), which.min(strata.mean),
        which(strata.max == min(strata.max[I])), which(strata.w.mean == min(strata.w.mean[I])))
  best.I = rbind(c(sigma[S[1]], r[S[1]]), c(sigma[S[2]], r[S[2]]),
                 c(sigma[S[3]], r[S[3]]), c(sigma[S[4]], r[S[4]]))
  colnames(best.I) = c("sigma","r")
  best = list(Fitted.strata.diff(X1, Z,
                                 cbps.k$cbps.fit[[best.I[1,1]]][[best.I[1,2]]]$propensity.score[,1],
                                 cbps.k$cbps.fit[[best.I[1,1]]][[best.I[1,2]]]$propensity.score[,2], ej),
              Fitted.strata.diff(X1, Z,
                                 cbps.k$cbps.fit[[best.I[2,1]]][[best.I[2,2]]]$propensity.score[,1],
                                 cbps.k$cbps.fit[[best.I[2,1]]][[best.I[2,2]]]$propensity.score[,2], ej),
              Fitted.strata.diff(X1, Z,
                                 cbps.k$cbps.fit[[best.I[3,1]]][[best.I[3,2]]]$propensity.score[,1],
                                 cbps.k$cbps.fit[[best.I[3,1]]][[best.I[3,2]]]$propensity.score[,2], ej),
              Fitted.strata.diff(X1, Z,
                                 cbps.k$cbps.fit[[best.I[4,1]]][[best.I[4,2]]]$propensity.score[,1],
                                 cbps.k$cbps.fit[[best.I[4,1]]][[best.I[4,2]]]$propensity.score[,2], ej))
  best.ps = cbind(cbps.k$cbps.fit[[best.I[1,1]]][[best.I[1,2]]]$propensity.score[,2],
                  cbps.k$cbps.fit[[best.I[2,1]]][[best.I[2,2]]]$propensity.score[,2],
                  cbps.k$cbps.fit[[best.I[3,1]]][[best.I[3,2]]]$propensity.score[,2],
                  cbps.k$cbps.fit[[best.I[4,1]]][[best.I[4,2]]]$propensity.score[,2],
                  cbps.k$cbps.fit[[best.I[1,1]]][[best.I[1,2]]]$propensity.score[,1],
                  cbps.k$cbps.fit[[best.I[2,1]]][[best.I[2,2]]]$propensity.score[,1],
                  cbps.k$cbps.fit[[best.I[3,1]]][[best.I[3,2]]]$propensity.score[,1],
                  cbps.k$cbps.fit[[best.I[4,1]]][[best.I[4,2]]]$propensity.score[,1])
  colnames(best.ps) = c("SO.mean.w","strata.mean.w","strata.max.w","strata.weight.mean.w",
                        "SO.ps","strata.mean.ps","strata.max.ps","strata.weight.mean.ps")
  out = list(best.index = best.I, best = best, best.ps = best.ps, all.max = all.max, strata.mean = strata.mean,
             strata.max = strata.max, strata.w.mean = strata.w.mean)
  return(out)
}



