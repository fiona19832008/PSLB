#' Kernel PCA using Gaussian kernel
#'
#' This function transforms the imput covariate space to the feature space using kernel PCA.
#' The kernel matrix K is constructed by Gaussian kernel. The input is the original covariate matrix and the number of candidate sigma from the kernel function used.
#' The function returns the list of the feature space and other parameters from the kernel PCA.
#'
#' @param X The covariate matrix with the rows corresponding to the subjects/participants and the columns corresponding to the covariates.
#' @param n_sigma The number of bandwidth value of the Gaussian kernel function used to generate the feature space. See details below. Default to 20.
#' @param k The number of top K eigen value calculated. Note that k should be smaller than the sample size. The minimum of k and the sample size will be used. Defaults to 500.
#' @details Kernel PCA is a simple way of implemention of RKHS to extend certain parametric methods with flexible modeling. Kernel PCA can extend the input covariates to a feature space of non-linear transformations defined by a kernel function.
#'          The kernel function (Gaussian kernel is implemented here) defines an inner product in the feature space between each pair of subjects.
#'          The kernel matrix K is constructed by calculating a measure of similarity between any two subjects with the Gaussian kernel function. And different sigma, the bandwitdh of the Gaussian kernel, determines the set of functions in the transformed space.
#'          Then eigen decomposition is conducted on the kernel matrix. The eigen vector matrix is the transformed feature matrix. For computational efficiency, only the top k eigen values are calculated.
#'          The range of the bandwidth sigma is the 0.1 and 0.9 quantiles of the L_2 distances between samples.
#'          The set of sigma includes n_sigma values from the function input within this range that are equally spaced on the log-scale.
#'          For each sigma value, the number of the eigen vectors corresponding to 90%, 95%, and 99% variance explained by the eigen values are calculated.
#' @return A list containing the following components:
#'         \itemize{
#'           \item{"features"}: A list of transformed feature space corresponding to each bandwidth sigma.
#'            Each elements of this list is a list whose elements are the feature space matrixes
#'            containing the number of the selected columns correspond to 90%, 95%, and 99% variance explained by the eigen values.
#'           \item{"r"}: A matrix containing the number of the eigen vectors corresponding to 90%, 95%, and 99% variance explained by the eigen values for each bandwith sigma.
#'           \item{"d"}: A list of eigen value corresponding to each bandwidth sigma.
#'           \item{"sigma_seq"}: The value of the bandwith sigma from the Gaussian kernel function.
#'         }
#'
#' @references Chen, W. and Zhang, H. (2007) The condition of kernelizing an algorithm and an equivalence between
#'             kernel methods. In Iberian Conference on Pattern Recognition and Image Analysis, 228-245. Springer.
#'
#' @examples
#' KS = Kang_Schafer_Simulation(n = 1000, seeds = 5050)
#' # Misspecified propensity score model
#' X = KS$Data[,7:10]
#' # Transform X with Gaussian kernel using 10 sigma
#' x.k = Gaussian_Kernel_feature(X = X, n_sigma = 10)
#'
#' @export
Gaussian_Kernel_feature = function(X, n_sigma = 20, k = 500) {
  n = nrow(X)
  k = min(k, n)
  kernel.fit = Gaussian_Kernel_choose_sigma(X, n_sigma = n_sigma, k = k)
  r99 = c()
  for(i in 1:length(kernel.fit$features)){
    d_cum = cumsum(kernel.fit$d[[i]])/sum(kernel.fit$d[[i]])
    r99[i] = sum(d_cum<=0.99)
  }
  r.m = max(r99, na.rm = TRUE)
  if (r.m <= k && abs(r.m - k) > 50) {
    kernel.fit = kernel.fit
  } else if (k == n) {
    kernel.fit = kernel.fit
  } else if (k <= 200) {
    kernel.fit = kernel.fit
  } else {
    kernel.fit = Gaussian_Kernel_choose_sigma(X, n_sigma = n_sigma, k = (r.m + 100))
  }
  features.out = list()
  r.out = matrix(nrow = length(kernel.fit$features), ncol = 3)
  for(i in 1:length(kernel.fit$features)) {
    d_cum = cumsum(kernel.fit$d[[i]])/sum(kernel.fit$d[[i]])
    r.med2 = sum(d_cum<=0.90)
    r.med3 = sum(d_cum<=0.95)
    r.max = sum(d_cum<=0.99)
    r = c(r.med2, r.med3, r.max)
    print(paste0("r=", r))
    r.out[i,] = r
    features.out[[i]] = list()
    for(j in 1:length(r)) {
      features.out[[i]][[j]] = kernel.fit$features[[i]][,1:r[j]]
    }
  }
  rownames(r.out) = paste0("sigma",seq(1:n_sigma))
  colnames(r.out) = c("0.90","0.95","0.99")
  out = list(features = features.out, r = r.out,
             d = kernel.fit$d, sigma_seq = kernel.fit$sigma_seq)
  return(out)
}

Gaussian_Kernel_choose_sigma = function(X, n_sigma, k = 500) {
  require(kernlab)
  require(RSpectra)
  set.seed(5050)
  X_mean = colMeans(X)
  X_sd = apply(X, 2, sd)
  X_std = matrix(nrow = nrow(X), ncol = ncol(X), byrow = T)
  for (i in 1:nrow(X)) {
    X_std[i, ] = (X[i, ] - X_mean)/X_sd
  }
  sigma = sigest(X_std)
  sigma_range = c(as.numeric(sigma[1]),  as.numeric(sigma[3]))
  log_sigma_range = log(sigma_range)
  sigma_seq = seq(from = log_sigma_range[1], to = log_sigma_range[2], length.out =  n_sigma)
  sigma_seq = exp(sigma_seq)
  eigen.out = list()
  d_out = list()
  features_out = list()
  for(i in 1:length(sigma_seq)) {
    rbf = rbfdot(sigma = sigma_seq[i])
    K = kernelMatrix(rbf, X_std)
    #K.eigen = eigen(K)
    #eigen = K.eigen$vectors
    defaultW <- getOption("warn")
    options(warn = -1)
    K.eigen = eigs(K, k = k)
    options(warn = defaultW)
    eigen = K.eigen$vectors
    d = K.eigen$values
    I = which(abs(d) <= 1e-10)
    d[I] = 0
    features = t(t(eigen)*sqrt(d))
    #features = features[, 1:r]
    #d = 1/d[1:r]
    eigen.out[[i]] = eigen
    d_out[[i]] = d
    #d_cum = cumsum(d)
    #d_cum_out[[i]] = d_cum/sum(d)
    features_out[[i]] = features
  }
  return(list(features = features_out, eigen.vector = eigen.out, d = d_out, sigma_seq = sigma_seq))
}


Form.formula = function(N) {
  out = paste0("Z~")
  if (length(N) == 1) {
    out = paste0(out, N[1])
  } else if (length(N) > 1) {
    out = paste0(out, N[1])
    for(i in 2:length(N)) {
      out = paste0(out, "+", N[i])
    }
  }
  return(out)
}

#' Kang and Schafer simulation
#'
#' This function generate the simulation scenarios presented in Kang and Schafer (2007)
#'
#' @param n The number of sample size.
#' @param beta The coefficient of the true propensity score model. Default to c(-1,0.5,-0.25,-0.1).
#' @param alpha The coefficient of the true outcome model. Defaults to c(210,27.4,13.7,13.7,13.7).
#' @param mu The mean of the covariates Z in the true propensity score model. Default to rep(0,4)
#' @param sd The variance of the covariates Z in the true propensity score model. Default to diag(4)
#' @param seeds The seed number
#' @return A list containing the following components:
#'         \itemize{
#'           \item{"Data"}: The simulated data matrix includes the outcome Y (1st column), the treatment assignment Tr (2nd column),
#'                          the covariates Z in the true propensit score mode (3rd to 6th column), the observed covariates X (7th to 10th column),
#'                          and the true propensity score PS (11th column).
#'           \item{"Treat.effect"}: The mean of the outcome Y
#'         }
#'
#' @references Kang, J. D. and Shafer, J. L. (2007) Demystifying double robustness: A comparison of alternative strategies for estimating a population mean from incomplete data. Statistical Science, 22, 523-539.
#' @export
Kang_Schafer_Simulation = function(n, beta = c(-1,0.5,-0.25,-0.1),
                                   alpha = c(210,27.4,13.7,13.7,13.7),
                                   mu = rep(0,4), sd = diag(4), seeds){
  require(MASS)
  set.seed(seeds)
  beta = matrix(beta, ncol = 1)
  alpha = matrix(alpha, ncol = 1)
  Z = mvrnorm(n, mu = mu, Sigma = sd)
  Z1 = cbind(rep(1,n),Z)
  epsilon = rnorm(n, mean = 0, sd = 1)
  Y = Z1%*%alpha + epsilon
  ps = exp(Z%*%beta)/(1+exp(Z%*%beta))
  D = c()
  for(i in 1:n) {
    D[i] = rbinom(1, size = 1, prob = ps[i])
  }
  X1 = exp(Z[,1]/2)
  X2 = Z[,2]/(1+exp(Z[,1])) + 10
  X3 = (((Z[,1]*Z[,3])/25)+0.6)^3
  X4 = (Z[,2] + Z[,4] + 20)^2
  X = cbind(X1, X2, X3, X4)
  out = cbind(Y, D, Z, X, ps)
  colnames(out) = c("Y", "Tr", "Z1", "Z2", "Z3", "Z4", "X1", "X2", "X3", "X4","PS")
  #treatment = mean(out[which(D==1),1]) - mean(out[which(D==0),1])
  treatment = mean(out[,1])
  return(list(Data = out, Treat.effect = treatment))
}


