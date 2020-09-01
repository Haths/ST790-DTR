#------------------------------------------------------------------------------#
# stratification estimator of tx effect
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
#------------------------------------------------------------------------------#
delta_S <- function(moPS, data, y, K) {

  #### Propensity Score Regression

  fitPS <- modelObj::fit(object = moPS, data = data, response = data$A)

  # estimated propensity score
  p1 <- modelObj::predict(object = fitPS, newdata = data) 
  if (is.matrix(x = p1)) p1 <- p1[,ncol(x = p1), drop = TRUE]

  #### Identify Strata

  # cutoff points for K groups 
  probs <- seq(from = 0.0, to = 1.0, length.out = K+1L)
  coff <- stats::quantile(x = p1, probs = probs) 
  coff[1L] <- 0.0
  coff[K+1L] <- 1.0

  # group ids for each individual
  grp <- colSums(x = sapply(X = p1, FUN = function(x){x <= coff}))

  #### Treatment Effect

  EY <- matrix(data = 0.0, nrow = K, ncol = 2L, 
               dimnames = list(NULL, c("0","1")))

  delta <- 0.0

  deltaj <- array(data = 0.0, dim = K, dimnames = list(1L:K))

  for (j in 1L:K) {

    gji <- grp == j
    
    EY[j,2L] <- mean(x = data$A[gji] * y[gji]) / mean(x = data$A[gji])
    EY[j,1L] <- mean({1.0-data$A[gji]} * y[gji]) / mean(x = 1.0 - data$A[gji])

    deltaj[j] <- EY[j,2L] - EY[j,1L]

    delta <- delta + sum(gji) / nrow(x = data) * deltaj[j]
  }

  delta <- unname(obj = delta)

  return( list(    "fitPS" = fitPS,
                      "cj" = coff,
                     "grp" = table(bin = grp),
                      "EY" = EY,
               "deltaHatj" = deltaj,
                "deltaHat" = delta) )

}
