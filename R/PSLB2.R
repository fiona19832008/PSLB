#' Propensity Score Analysis with Local Balance (PSLB) Estimation
#' 
#' PSLB nonparametrically estimates the propensity score weights that achieve covariate balance in both whole study population (global balance) 
#' and the propensity score stratified (PS-stratified) sub-populations (local balance). Therefore, the PSLB estimated propensity score are appropriate 
#' since they approximates the balancing score, when the propensity score model is subject to model misspecification. PSLB include two steps (PSLB 1 and PSLB 2).
#' In the 1st step, PSLB 1 implements a flexible form of covariate balancing propensity score (CBPS), and finds the model with the best local balance 
#' among the model pool. In the 2nd step, PSLB 2 refines the score from PSLB 1 so that the estimated score approximately balances the covariates to a 
#' pre-specified level in the local neighborhoods (each neighborhood, a small inverval between 0 and 1, represents a sub-population whose estimated score falls in this inverval). 
#' Note that PSLB 2 can be used to improve the local balance of propensity score estimated by any method, such as CBPS (see the examples). 
#' 
#' @seealso [PSLB1()] for the detailed description of PSLB 1.
#' 
#' @param X The covariate matrix with the rows corresponding to the subjects/participants and the columns corresponding to the covariates.
#'          The covariate matrix X can not contain missing value. The function will stop if NA value is detected in X.
#' @param Z The binary treatment indicator. A vector with 2 unique numeric values in 0 = untreated and 1 = treated.
#'          Its length should be the same as the number of all subjects/participants without missing value.
#' @param n_sigma The number of bandwidth value of the Gaussian kernel function used to generate the feature space in PSLB 1.
#'                See details in the description of function \code{\link{Gaussian_Kernel_feature}} and \code{\link{PSLB1}}. Default to 20.
#' @param ej The local neighborhoods for PSLB 1 and PSLB 2. It can be one matrix or a list of two matrixes containing values greater or equal to 0 and less or equal to 1.
#'           If ej is one matrix, PSLB 1 and PSLB 2 both use this matrix as their local neighborhoods. 
#'           If ej is a list of matrixes with length of two, PSLB 1 uses the first element of the list and PSLB 2 uses the second
#'           element of the list as their local neighborhoods. The rows of the matrixes represent the local neighborhoods. 
#'           There are two columns, where the first columns corresponding the start and the second column corresponding the end point of the local neighborhoods. 
#' @param selectInX The matrix used to evaluate the local balance when tune the model in PSLB 1. selectInX = "cov" by default, which evaluates the local covariate balance on the
#'                  input covariate matrix X. If selectInX = "feature", the local covariate balance is evaluated on the feature space transformed by kernel PCA.
#' @param k The number of top K eigen value calculated in PSLB 1. Note that k should be smaller than the sample size. The minimum of k and the sample size will be used.
#'          Defaults to 500. See details in the description of function \code{\link{Gaussian_Kernel_feature}} and \code{\link{PSLB1}}.
#' @param method Choose "exact" to fit the justidentified CBPS model in PSLB 1; choose "over" to fit the overidentified CBPS model in PSLB 1. 
#'               See details in the description of function \code{\link{PSLB1}}. Default is "exact".
#' @param standardize Default is TRUE, which normalizes weights of PSLB 1 to sum to 1 within treated/untreated group.
#'                    Set to FALSE to return IPW weights for ATE.
#' @param criteria The statistics used in PSLB 1 to select the model with the optimized local balance. 
#'                 Choose "mean" to use the mean of the absolute S/D among all covariates and neighborhoods for the tuning parameters by T4LB;
#'                 choose "max" to use the max of the absolute S/D among all covariates and neighborhoods for the tuning parameters by T4LB.
#'                 Default is "mean".
#' @param p.hat User supplied propensity score. When p.hat is provided, PSLB 1 will not be fitted. PSLB 2 uses the provided p.hat to determine the 
#'              PS-stratified sub-populations. Then IPW weight is estimated such that the propensity score within each strata are adjusted to 
#'              minimizing the local imbalance in that interval while holding the PS-stratified sub-populations unchanged from p.hat.
#'              When p.hat = NULL, the function fits the PSLB 1 and PSLB 2 subsequently as discussed in the description. The default is NULL.
#' @param epsilon A list of factors for controlling the local covariate imbalance parameters
#' \itemize{
#'   \item{"input"}: The matrix of candidate local imbalance parameters of choice. If input = NULL, the imbalance parameters can be set 
#'                   by our default rules. However, the users can supply their own values. The "input" must be a matrix with rows corresponding to
#'                   each set of candidate imbalance parameters and columns corresponding to the local neighborhoods. Therefor, its number of columns
#'                   must equal to the number of local neighborhoods.
#'   \item{"e.min"}: The minimum candidate value of the imbalance parameters in log2 scale. Default is -10
#'   \item{"bylength"}: The step size of the candidate imbalance parameters in log2 scale. Default is 0.2
#'   \item{"a"}: The factor determines the maximum candidate value of the imbalance parameters. The largest candidate value is set to be
#'               "a" times of the estimated standard deviation of the weighted mean difference of covariates between the two treatment groups in
#'               each sub-population. The default value of a is 3. When PSLB 2 is unsolvable, increasing "a" may result in a solution. However,
#'               this solution gives larger local balance compared to when smaller "a" is used.
#'   \item{"epsilon.div"}: The vector groups the neighboring stratas who share the same value of epsilon. When use non-overlapping strata (local neighborhood), 
#'                         "epsilon.div" is unnecessary. When use overlapping strata (local neighborhood), "epsilon.div" must be provided. Each element of "epsilon.div" 
#'                         corresponds the number of neighboring strata. For example, epsilon.div = c(2,3,2) means that the same imbalance parameter (epsilon) is used in 
#'                         strata 1 and 2; strata 3, 4, and 5; or strata 6 and 7, respectively. If the user would like to have different imbalance parameter in each neighborhood,
#'                         "epsilon.div" is a vector of 1s, whose length equals to the number of stratas. For example, if there is 5 local neighborhoods, epsilon.div = rep(1,5)
#'                         will allow different epsilon value in each strata.
#'                }
#' 
#' @return A list containing the following components:
#'         \itemize{ 
#'           \item{"PSLB2"}: A list containing the following components:
#'               \itemize{ 
#'                 \item{"weight"}: The optimal weights estimated by PSLB 2. If the fitting of PSLB 2 is unsuccess, weight returns NULL value.
#'                 \item{"local.sample.size"}: The matrix contains the sample size in total (n), in treated group (n1), and in untreated group (n0) for each local
#'                                             neighborhoods. Each row of this matrix represents each neighborhood. The columns corresponds to n, n1/n, n0 and n1.
#'                 \item{"p.hat"}: The propensity score used to determine the PS-stratified sub-population in PSLB 2. If PSLB 1 is fitted, the output "p.hat" is the
#'                                 estimated propensity score from PSLB 1. If PSLB 1 is not fitted, the output "p.hat" is just the input "p.hat".
#'                 \item{"balance"}: The matrix containing the global and local balance of the estimated weights from PSLB 2. The first row contains the absolute S/D of each covariates
#'                                   in the whole study population, which represents the global covariate balance. The rest rows contain the absolute
#'                                   S/D of each covariates in the PS-stratified sub-population corresponding to the loal neighborhoods of PSLB 2, which
#'                                   represent the local covariate balance. If the fitting of PSLB 2 is unsuccess, balance returns NULL value.
#'                 \item{"epsilon"}: The selected local imbalance parameter for each local neighborhood.
#'                 \item{"epsilon.all"}: All the candidate values of the local imbalance parameter epsilon. If use non-overlapping strata, "epsilon.all" is a list whose
#'                                       each element corresponds to each strata. If use overlapping strata, "epsilon.all" is a vector since all local neighborhoods share
#'                                       one set of candidate values of epsilon. 
#'               }
#'           \item{"PSLB1"}: The output list of function \code{\link{PSLB1}}
#'           \item{"status"}: An integer code. 0 indicates successful completion of PSLB 2 estimation. Possible error codes are
#'                            1 and 2. 1 indicates that PSLB 2 can not find a solution under the current covariate imbalance condition.
#'                            2 indicates that there is an error occured when fits the PSLB2.
#'           \item{"message"}: A character string giving the fitting status of PSLB 2.
#'           \item{"X"}: The input covariate matrix
#'           \item{"treat"}: The treatment assignment vector used
#'         }
#' 
#' @references Li, Y. and Li, L. Propensity score analysis with local balance. manuscript under preperation. 
#' 
#' @examples 
#' KS = Kang_Schafer_Simulation(n = 1000, seeds = 371834)
#' # Misspecified propensity score model
#' X = KS$Data[,7:10]
#' Z = KS$Data[,2]
#' # Specify the local neighborhoods
#' ej1 = cbind(seq(0,0.8,0.2),seq(0.2,1,0.2)) # non-overlapping strata
#' ej2 = cbind(seq(0,0.6,0.1),seq(0.4,1,0.1)) # overlapping strata
#' ej = list(ej1, ej2)
#' print(ej1)
#' print(ej2)
#' # Fit the PSLB model with overlapping strata
#' PSLB.fit = PSLB(X = X, Z = Z, n_sigma = 10, ej = ej, p.hat = NULL, epsilon = list(epsilon.div = c(2,3,2)))
#' str(PSLB.fit)
#' # Global and local balance of PSLB 2 in ej2
#' round(PSLB.fit$PSLB2$balance,3)
#' # Global and local balance of PSLB 1 in ej2
#' round(abs(Fitted.strata.diff(X = X, Z = Z, ps = PSLB.fit$PSLB1$propensity.score, ej = ej2)),3)
#' # Fit the PSLB model with non-overlapping strata
#' PSLB.fit.ej1 = PSLB(X = X, Z = Z, n_sigma = 10, ej = ej1, p.hat = PSLB.fit$PSLB1$propensity.score, epsilon = list(a = 2))
#' # Global and local balance of PSLB 2 in ej1
#' round(PSLB.fit.ej1$PSLB2$balance, 3)
#' # Fit the CBPS model
#' cbps.fit = CBPS(Z~X, ATT = 0, method = "exact")
#' cbps.ps = cbps.fit$fitted.values
#' # Improve the local balance of the estimated propensity score from CBPS
#' cbps.pslb.fit = PSLB(X = X, Z = Z, ej = ej, p.hat = cbps.ps, epsilon = list(epsilon.div = c(2,3,2)))
#' round(cbps.pslb.fit$PSLB2$balance,3)
#' # Global and local imbalance of CBPS
#' cbps.local = abs(Fitted.strata.diff(X = X, Z = Z, ps = cbps.ps, ej = ej2))
#' round(cbps.local, 3)
#' 
#' @export
PSLB = function(X, Z, n_sigma = 20, ej, selectInX = "cov", k = 500, method = "exact", standardize = TRUE,
                criteria = "mean", p.hat = NULL, epsilon = list(input = NULL, e.min = -10, bylength = 0.2, a = 3, epsilon.div = NULL)) {
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
  # Local neighborhoods for PSLB 1 (ej1) and PSLB 2 (ej2)
  if (is.list(ej)) {
    ej1 = ej[[1]]
    ej2 = ej[[2]]
  } else if (is.matrix(ej)) {
    ej1 = ej
    ej2 = ej
  }
  # Find if the strata for PSLB 2 is overlapped or not
  over.I = OverlapOrNot(ej2)
  # over.I == 0: nonoverlapping strata; over.I != 0: overlapping strata
  if (over.I != 0 && is.null(epsilon$epsilon.div)) {
    stop("The vector that groups the local neighborhoods must be provided.")
  }
  # Assign values in epsilon list if not provided
  if (is.null(epsilon$e.min)) {
    epsilon$e.min = -10
  }
  if (is.null(epsilon$bylength)) {
    epsilon$bylength = 0.2
  }
  if (is.null(epsilon$a)) {
    epsilon$a = 3
  }
  # Step 1: Fit PSLB 1 and tune with T4LB
  if (is.null(p.hat[1])) {
    PSLB1.fit = PSLB1(X = X, Z = Z, n_sigma = n_sigma, ej = ej1, selectInX = selectInX, k = k,
                      method = method, standardize = standardize, criteria = criteria)
    if (is.null(PSLB1.fit[1])) {
      stop("The 1st step of PSLB algorithm (PSLB 1) failed.")
    } else {
      p.hat = PSLB1.fit$propensity.score
    }
  } else if (!is.null(p.hat[1])) {
    p.hat = p.hat
    PSLB1.fit = list(NA, "message" = "Did not fit PSLB 1.")
  }
  # Step 2: Fit PSLB 2 with non-overlapping or overlapping strata
  # non-overlapping strata: over.I == 0; overlapping strata: over.I!= 0
  if (over.I == 0) {
    PSLB2.fit = NULL
    tryCatch({PSLB2.fit = PSWquadra.fit(p.hat = p.hat, X = X, Z = Z, 
                                        epsilon = list(input = epsilon$input, e.min = epsilon$e.min, 
                                                       bylength = epsilon$bylength, a = epsilon$a),
                                        ej = ej2, sol = list(sol_name = "osqp", min.sub = 0.1))},
             error = function(e) {
               print("PSLB2 failed due to fitting error.")
             })
  } else if (over.I != 0) {
    # X_global is any input variable that only adjust global balance. 
    PSLB2.fit = NULL
    tryCatch({PSLB2.fit = PSWquadra.fit.OverlapStrata(p.hat = p.hat, X = X, Z = Z, 
                                        epsilon = list(input = epsilon$input, e.min = epsilon$e.min, 
                                                       bylength = epsilon$bylength, a = epsilon$a),
                                        ej = ej2, X_global = NULL, sol = list(sol_name = "osqp", epsilon.div = epsilon$epsilon.div, min.sub = 0.1))},
             error = function(e) {
               print("An error occured when fits the PSLB2")
             })
  }
  # If PSLB2 failed at fitting
  if (is.null(PSLB2.fit[1])) {
    # status 2: An error occured when fits the PSLB2; Please check if the inputs are appropriate.
    LSS = sample.size.calculate(Z = Z, p.hat = p.hat, ej = ej2)
    rownames(LSS) = paste0(ej2[,1], "-", ej2[,2])
    PSLB2 = list(weight = NULL, local.sample.size = LSS, p.hat = p.hat, balance = NULL, 
                 epsilon = NULL, epsilon.all = NULL)
    out = list(PSLB2 = PSLB2, PSLB1 = PSLB1.fit, status = 2,
               message = "An error occured when fits the PSLB2",
               X = X, treat = Z)
  } else if (!is.null(PSLB2.fit[1])) {
    if (sum(is.na(PSLB2.fit$weight)) != 0) {
      weight = NULL
      balance = NULL
      status = 1
      mg = "PSLB 2 can not find a solution under the current covariate imbalance condition."
    } else if (sum(is.na(PSLB2.fit$weight)) == 0) {
      weight = PSLB2.fit$weight
      balance = abs(Fitted.strata.diff(X = X, Z = Z, ps = p.hat, weight = PSLB2.fit$weight,
                                       ej = ej2))
      status = 0
      mg = "PSLB 2 fitted successfully."
    }
    e.out = PSLB2.fit$epsilon[,1:3]
    rownames(e.out) = paste0("strata",seq(1,nrow(ej2)))
    colnames(e.out) = c("start","end","epsilon")
    LSS = PSLB2.fit$sample.size
    rownames(LSS) = paste0(e.out[,1], "-", e.out[,2])
    PSLB2 = list(weight = weight, local.sample.size = LSS, p.hat = p.hat, balance = balance,
                 epsilon = e.out, epsilon.all = PSLB2.fit$epsilon.list)
    out = list(PSLB2 = PSLB2, PSLB1 = PSLB1.fit, status = status,
               message = mg, X = X, treat = Z)
  }
  return(out)
}

OverlapOrNot = function(ej) {
  p = nrow(ej)
  start = ej[-1,1]
  end = ej[-p,2]
  I = start < end
  out = sum(I)
  return(out)
}

