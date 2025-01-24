% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PSWweight.R
\name{PS.ConvertFrom.Weight}
\alias{PS.ConvertFrom.Weight}
\title{Convert the Inverse Probability Weight (IPW) for ATE to Propensity Score}
\usage{
PS.ConvertFrom.Weight(w, Z)
}
\arguments{
\item{w}{The input IPW weights.}

\item{Z}{The binary treatment indicator. A vector with 2 unique numeric values in 0 = untreated and 1 = treated.}
}
\value{
The vector contanining the converted propensity score.
}
\description{
This function converts the IPW weights from PSLB 2 to the corresponding propensity scores. When the subject i is
in the untreated group, the converted propensity score is 1-(1/weight_i). When the subject i is in the treated group,
converted propensity score is 1/weight_i.
}
