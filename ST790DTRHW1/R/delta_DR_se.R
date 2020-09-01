#------------------------------------------------------------------------------#
# doubly robust estimator of tx effect and its standard error
#
# REQUIRES
#   swv_ML(), swv_OLS(), and delta_DR_swv()
#
# ASSUMPTIONS
#   tx is denoted by A and is binary coded as {0,1}
#   outcome regression model is a linear model and parameters are estimated 
#     using OLS
#   propensity score regression model is a logistic model and parameters are
#     estimated using ML
#
# INPUTS
#  moOR : a modeling object specifying the outcome regression step
#         *** must be a linear model ***
#  moPS : a modeling object specifying the propensity score regression step
#         *** must be a logistic model ***
#  data : a data.frame containing covariates and tx
#         *** tx must be named 'A' and coded 0/1 ***
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
#  sigmaHat : the estimated standard error
#------------------------------------------------------------------------------#
delta_DR_se <- function(moOR, moPS, data, y) {

  #### ML components
  ML <- swv_ML(mo = moPS, data = data, y = data$A) 

  #### OLS components
  OLS <- swv_OLS(mo = moOR, data = data, y = y) 

  #### estimator components
  DR <- delta_DR_swv(moOR = moOR, moPS = moPS, data = data, y = y)

  #### 1,1 Component of Sandwich Estimator

  # ML contribution
  temp <- DR$dMdG %*% ML$invdM
  sig11ML <- temp %*% ML$MM %*% t(x = temp)

  # OLS contribution
  temp <- DR$dMdB %*% OLS$invdM
  sig11OLS <- temp %*% OLS$MM %*% t(x = temp)

  sig11 <- drop(x = DR$MM + sig11ML + sig11OLS)

  return( c(DR$delta, "sigmaHat" = sqrt(x = sig11 / nrow(x = data))) )

}
