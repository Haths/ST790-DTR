#------------------------------------------------------------------------------#
# outcome regression estimator of tx effect
#
# ASSUMPTIONS
#   tx is denoted by A and is binary coded as {0,1}
#
# INPUTS
#  moOR : a modeling object specifying the outcome regression step
#  data : a data.frame containing covariates and tx
#     y : a vector containing the outcome of interest
#
# RETURNS
#  a list containing
#     fitOR : a modelObjFit object containing the results of the
#             outcome regression step
#        EY : the sample mean of the outcome under each tx option
#  deltaHat : the outcome regression estimator for the tx effect
#------------------------------------------------------------------------------#
delta_OR <- function(moOR, data, y) {

  #### Outcome Regression 

  fitOR <- modelObj::fit(object = moOR, data = data, response = y)

  # Q(X,0;betaHat)
  data$A <- 0L 
  Q0 <- drop(x = modelObj::predict(object = fitOR, newdata = data))

  # Q(X,1;betaHat)
  data$A <- 1L 
  Q1 <- drop(x = modelObj::predict(object = fitOR, newdata = data))

  #### Tx Effect

  EY <- array(data = 0.0, dim = 2L, dimnames = list(c("0","1")))
  EY[2L] <- mean(x = Q1)
  EY[1L] <- mean(x = Q0)
  delta <- unname(obj = EY[2L] - EY[1L])

  return( list(   "fitOR" = fitOR, 
                     "EY" = EY,
               "deltaHat" = delta) )
}
