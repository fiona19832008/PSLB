% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gaussian.kernel.R
\name{Kang_Schafer_Simulation}
\alias{Kang_Schafer_Simulation}
\title{Kang and Schafer simulation}
\usage{
Kang_Schafer_Simulation(
  n,
  beta = c(-1, 0.5, -0.25, -0.1),
  alpha = c(210, 27.4, 13.7, 13.7, 13.7),
  mu = rep(0, 4),
  sd = diag(4),
  seeds
)
}
\arguments{
\item{n}{The number of sample size.}

\item{beta}{The coefficient of the true propensity score model. Default to c(-1,0.5,-0.25,-0.1).}

\item{alpha}{The coefficient of the true outcome model. Defaults to c(210,27.4,13.7,13.7,13.7).}

\item{mu}{The mean of the covariates Z in the true propensity score model. Default to rep(0,4)}

\item{sd}{The variance of the covariates Z in the true propensity score model. Default to diag(4)}

\item{seeds}{The seed number}
}
\value{
A list containing the following components:
\itemize{
\item{"Data"}: The simulated data matrix includes the outcome Y (1st column), the treatment assignment Tr (2nd column),
the covariates Z in the true propensit score mode (3rd to 6th column), the observed covariates X (7th to 10th column),
and the true propensity score PS (11th column).
\item{"Treat.effect"}: The mean of the outcome Y
}
}
\description{
This function generate the simulation scenarios presented in Kang and Schafer (2007)
}
\references{
Kang, J. D. and Shafer, J. L. (2007) Demystifying double robustness: A comparison of alternative strategies for estimating a population mean from incomplete data. Statistical Science, 22, 523-539.
}
