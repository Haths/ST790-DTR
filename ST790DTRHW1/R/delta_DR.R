#------------------------------------------------------------------------------#
# doubly robust estimator of tx effect 
#
# ASSUMPTIONS
#   tx is denoted by A and is binary coded as {0,1}
#
# INPUTS
#  moOR : a modeling object specifying the outcome regression step
#  moPS : a modeling object specifying the propensity score regression step
#  data : a data.frame containing covariates and tx
#     y : a vector containing the outcome of interest
#
# RETURNS
#  a list containing
#     fitOR : a modelObjFit object containing the results of the
#             outcome regression step
#     fitPS : a modelObjFit object containing the results of the
#             propensity score regression step
#        EY : the sample mean of the outcome under each tx option
#  deltaHat : the doubly robust estimator for the tx effect
#------------------------------------------------------------------------------#
delta_DR <- function(moOR, moPS, data, y) {

  #### Propensity Score

  fitPS <- modelObj::fit(object = moPS, data = data, response = data$A)

  # estimated propensity score
  p1 <- modelObj::predict(object = fitPS, newdata = data) 
  if (is.matrix(x = p1)) p1 <- p1[,ncol(x = p1), drop = TRUE]

  #### Outcome Regression 

  fitOR <- modelObj::fit(object = moOR, data = data, response = y)

  # store tx variable
  A <- data$A 

  # estimated Q-function when all A=0
  data$A <- 0L 
  Q0 <- drop(x = modelObj::predict(object = fitOR, newdata = data))

  # estimated Q-function when all A=1
  data$A <- 1L 
  Q1 <- drop(x = modelObj::predict(object = fitOR, newdata = data))

  #### Treatment Effect

  EY <- array(data = 0.0, dim = 2L, dimnames = list(c("0","1")))
  aug <- {A - p1} * {Q1 / p1 + Q0 / {1.0 - p1}}
  EY[2L] <- mean(x = {A == 1L} * {y / p1 - aug})
  EY[1L] <- mean(x = {A == 0L} * {y / {1.0 - p1} + aug})
  delta <- unname(obj = EY[2L] - EY[1L])

  return( list(   "fitOR" = fitOR,
                  "fitPS" = fitPS,
               "deltaHat" = delta,
                     "EY" = EY) )

}
