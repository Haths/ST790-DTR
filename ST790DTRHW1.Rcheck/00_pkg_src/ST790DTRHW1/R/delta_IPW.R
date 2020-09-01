#------------------------------------------------------------------------------#
# IPW estimator of tx effect 
#
# ASSUMPTIONS
#   tx is denoted by A and is binary coded as {0,1}
#
# INPUTS
#  moPS : a modeling object specifying the propensity score regression step
#  data : a data.frame containing covariates and tx
#     y : a vector containing the outcome of interest
#
# RETURNS
#  a list containing
#     fitPS : a modelObjFit object containing the results of the
#             propensity score regression step
#        EY : the sample mean of the outcome under each tx option
#  deltaHat : the outcome regression estimator for the tx effect
#------------------------------------------------------------------------------#
delta_IPW <- function(moPS, data, y) {

  #### Propensity Score

  fitPS <- modelObj::fit(object = moPS, data = data, response = data$A)

  # estimated propensity score
  p1 <- modelObj::predict(object = fitPS, newdata = data) 
  if (is.matrix(x = p1)) p1 <- p1[,ncol(x = p1), drop = TRUE]

  #### Treatment Effect
  
  EY <- array(data = 0.0, dim = 2L, dimnames = list(c("0","1")))
  EY[1L] <- mean(x = {1.0 - data$A} * y / {1.0 - p1})
  EY[2L] <- mean(x = data$A * y / p1) 
  delta <- unname(obj = EY[2L] - EY[1L])

  return( list(   "fitPS" = fitPS,
                     "EY" = EY,
               "deltaHat" = delta) )
}
