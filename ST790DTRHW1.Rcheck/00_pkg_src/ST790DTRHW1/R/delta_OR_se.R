#------------------------------------------------------------------------------#
# outcome regression estimator of tx effect and its standard error
#
# REQUIRES
#   swv_OLS() and delta_OR_swv()
#
# ASSUMPTIONS
#   tx is denoted by A and is binary coded as {0,1}
#   outcome regression model is a linear model and parameters are estimated 
#     using OLS
#
# INPUTS
#  moOR : a modeling object specifying the outcome regression step
#         *** must be a linear model ***
#  data : a data.frame containing covariates and tx
#         *** tx must be named 'A' and coded 0/1 ***
#     y : a vector containing the outcome of interest
#
# RETURNS
#  a list containing
#  deltaHat : the outcome regression estimator for the tx effect
#        EY : the sample mean of the outcome under each tx option
#     fitOR : a modelObjFit object containing the results of the
#             outcome regression step
#  sigmaHat : the estimated standard error
#------------------------------------------------------------------------------#
delta_OR_se <- function(moOR, data, y) {

  #### OLS components
  OLS <- swv_OLS(mo = moOR, data = data, y = y) 

  #### estimator components
  OR <- delta_OR_swv(moOR = moOR, data = data, y = y)

  #### 1,1 Component of Sandwich Estimator

  # OLS contribution
  temp <- OR$dMdB %*% OLS$invdM
  sig11OLS <- temp %*% OLS$MM %*% t(x = temp)

  sig11 <- drop(x = OR$MM + sig11OLS)

  return( c(OR$delta, "sigmaHat" = sqrt(x = sig11 / nrow(x = data))) )

}
