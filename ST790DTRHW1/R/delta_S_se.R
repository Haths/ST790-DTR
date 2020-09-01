#------------------------------------------------------------------------------#
# stratification estimator of tx effect and its standard error
#
# REQUIRES
#   delta_S()
#
# ASSUMPTIONS
#   tx is denoted by A and is binary coded as {0,1}
#
# INPUTS
#  moPS : a modeling object specifying the propensity score regression step
#  data : a data.frame containing covariates and tx
#         *** tx must be named 'A' and coded 0/1 ***
#     y : a vector containing the outcome of interest
#     K : the integer number of strata
#
# RETURNS
#  a list containing
#      fitPS : a modelObjFit object containing the results of the
#              propensity score regression step
#         cj : the propensity score boundaries defining the strata
#        grp : the number of individuals in each stratum
#         EY : the sample mean of the outcome under each tx option
#  deltaHatj : the stratification estimator for the tx effect for each strata
#   deltaHat : the stratification estimator for the tx effect
#   sigmaHat : the estimated standard error
#------------------------------------------------------------------------------#
delta_S_se <- function(moPS, data, y, K) {

  #### Treatment Effect

  delta <- delta_S(moPS = moPS, data = data, y = y, K = K)

  #### Standard Error

  # estimated propensity score
  p1 <- modelObj::predict(object = delta$fitPS, newdata = data) 
  if (is.matrix(x = p1)) p1 <- p1[,ncol(x = p1), drop = TRUE]

  # group ids for each individual
  grp <- colSums(x = sapply(X = p1, FUN = function(x){x <= delta$cj}))

  sigmaSq <- 0.0

  for (j in 1L:K) {

    gji <- grp == j
    
    Y1ji <- data$A[gji] * y[gji] / mean(x = data$A[gji])
    Y0ji <- {1.0-data$A[gji]} * y[gji] / mean(x = 1.0 - data$A[gji])

    deltaji <- {Y1ji - Y0ji}

    sigmaSq <- sigmaSq + mean(x = {deltaji - delta$deltaHatj[j]}^2)
  }

  sigmaSq <- sigmaSq / K

  return( c(delta, "sigmaHat" = sqrt(x = sigmaSq / nrow(x = data))) )
                
}
