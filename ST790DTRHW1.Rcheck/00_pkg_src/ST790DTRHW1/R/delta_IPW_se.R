#------------------------------------------------------------------------------#
# IPW estimator of the tx effect and its standard error
#
# REQUIRES
#   swv_ML() and delta_IPW_swv()
#
# ASSUMPTIONS
#   tx is denoted by A and is binary coded as {0,1}
#   propensity score regression model is a logistic model and parameters are
#     estimated using ML
#
# INPUTS
#  moOR : a modeling object specifying the propensity score regression step
#         *** must be a logistic model ***
#  data : a data.frame containing covariates and tx
#         *** tx must be named 'A' and coded 0/1 ***
#     y : a vector containing the outcome of interest
#
# RETURNS
#  a list containing
#     fitPS : a modelObjFit object containing the results of the
#             propensity score regression step
#        EY : the sample mean of the outcome under each tx option
#  deltaHat : the outcome regression estimator for the tx effect
#  sigmaHat : the estimated standard error
#------------------------------------------------------------------------------#
delta_IPW_se <- function(moPS, data, y){

  #### ML components
  ML <- swv_ML(mo = moPS, data = data, y = data$A) 

  #### estimator components
  IPW <- delta_IPW_swv(moPS = moPS, data = data, y = y)

  #### 1,1 Component of Sandwich Estimator

  # ML contribution
  temp <- IPW$dMdG %*% ML$invdM
  sig11ML <- temp %*% ML$MM %*% t(x = temp)

  sig11 <- drop(x = IPW$MM + sig11ML)

  return( c(IPW$delta, "sigmaHat" = sqrt(x = sig11 / nrow(x = data))) )

}
