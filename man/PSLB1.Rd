% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CBPS.kernel.R
\name{PSLB1}
\alias{PSLB1}
\title{Propensity Score Analysis with Local Balance 1 (PSLB 1) Estimation}
\usage{
PSLB1(
  X,
  Z,
  n_sigma = 20,
  ej,
  selectInX = "cov",
  k = 500,
  method = "exact",
  standardize = TRUE,
  criteria = "mean"
)
}
\arguments{
\item{X}{The covariate matrix with the rows corresponding to the subjects/participants and the columns corresponding to the covariates.
The covariate matrix X can not contain missing value. The function will stop if NA value is detected in X.}

\item{Z}{The binary treatment indicator. A vector with 2 unique numeric values in 0 = untreated and 1 = treated.
Its length should be the same as the number of all subjects/participants without missing value.
The function will stop if NA value is detected in Z.}

\item{n_sigma}{The number of bandwidth value of the Gaussian kernel function used to generate the feature space.
See details in the description of function \code{\link{Gaussian_Kernel_feature}}. Default to 20.}

\item{ej}{The matrix of the local neighborhoods with its rows representing the neighborhoods. It contains two columns of values greater or equal to 0 and less or equal to 1.
The first columns are the start and the second column are the end point of the local neighborhoods.}

\item{selectInX}{The matrix used to evaluate the local balance. selectInX = "cov" by default, which evaluates the local covariate balance on the
input covariate matrix X. If selectInX = "feature", the local covariate balance is evaluated on the feature space transformed by kernel PCA.}

\item{k}{The number of top K eigen value calculated. Note that k should be smaller than the sample size. The minimum of k and the sample size will be used.
Defaults to 500. See details in the description of function \code{\link{Gaussian_Kernel_feature}}.}

\item{method}{Choose "exact" to fit the justidentified CBPS model; choose "over" to fit the overidentified CBPS model. See details in Imai and Ratkovic (2014).
Default is "exact".}

\item{standardize}{Default is TRUE, which normalizes weights to sum to 1 within treated/untreated group.
Set to FALSE to return IPW weights for ATE.}

\item{criteria}{Choose "mean" to use the mean of the absolute S/D among all covariates and neighborhoods for the tuning parameters by T4LB;
choose "max" to use the max of the absolute S/D among all covariates and neighborhoods for the tuning parameters by T4LB.
Default is "mean".}
}
\value{
A list containing the following components:
\itemize{
\item{"weight"}: The IPW weights for binary ATE calculated by the fitted propensity score given by Tr/ps + (1-Tr)/(1-ps),
where Tr is the treatment assignment and ps is the fitted propensity score. This expression for weight
is before standardization (i.e. with standardize = FALSE). Standardization will make weights sum to 1
within untreated/treated group.
\item{"propensity.score"}: The fitted propensity score
\item{"balance"}: The matrix containing the global and local balance of the estimated propensity score. The first row contains the absolute S/D of each covariates
in the whole study population, which represents the global covariate balance. The rest rows contain the absolute
S/D of each covariates in the PS-stratified sub-population corresponding to the loal neighborhoods of "ej", which
represent the local covariate balance.
\item{"coefficients"}: The coefficient of the fitted propensity score model
\item{"feature"}: The selected features by kernel PCA with the best loal covariate balance.
\item{"best.para"}: The parameters (bandwidth for Gaussian kernel and number of PCs) chosen with the best local covariate balance.
\item{"treat"}: The treatment assignment vector used
}
}
\description{
PSLB1 estimates propensity scores by implementing a flexible form of the covariate balancing propensity score (CBPS) using kernel PCA,
and tunes parameters (the bandwidth of the Gaussian kernel and the number of PCs) by the tuning for local balance (T4LB) algorithm, which finds
the model with the best local balance (minimum absolute standardized difference (S/D) in the PS-stratified sub-populations) among the model pool. The method
searches for the propensity score model with the best local balance while controling the global balance. The estimation is considered as "fail" if the minimum
absolute S/D of input covariates in the whole population is more than 10\% (uncontrolled global balance), and the function output NA
value with an error message. The local balance is evaluated by a statistic (mean or max) of the absolute S/D in the PS-stratified sub-populations.
The method only takes binary treatment.
}
\examples{
KS = Kang_Schafer_Simulation(n = 1500, seeds = 5050)
# Misspecified propensity score model
X = KS$Data[,7:10]
Z = KS$Data[,2]
# Specify the local neighborhoods
ej = cbind(seq(0,0.8,0.2),seq(0.2,1,0.2))
print(ej)
# PSLB 1 fitting
fit = PSLB1(X = X, Z = Z, n_sigma = 10, ej = ej)
print(fit$balance)

}
\references{
Imai, K. and Ratkovic, M. (2014) Covariate balancing propensity score. Journal of the
Royal Statistical Society: Series B (Statistical Methodology), 76, 243-263.
}
