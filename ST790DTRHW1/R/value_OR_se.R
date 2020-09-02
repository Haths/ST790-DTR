#------------------------------------------------------------------------------#
# value and standard error for the outcome regression estimator of the value
# of a fixed regime using the sandwich estimator of variance
#
# REQUIRES
#   swv_OLS() and value_OR_swv()
#
# ASSUMPTIONS
#   tx is binary coded as {0,1}
#   moOR is a linear model with parameters to be estimated using OLS
#
# INPUTS
#    moOR : modeling object specifying outcome regression
#           *** must be a linear model ***
#    data : data.frame containing covariates and tx
#           *** tx must be coded 0/1 ***
#       y : vector of outcome of interest
#  txName : tx column in data (tx must be coded as 0/1)
#  regime : 0/1 vector containing recommended tx
#
# RETURNS
#  a list containing
#     fitOR : value object returned by outcome regression
#        EY : mean response for each tx
#  valueHat : estimated value
#  sigmaHat : estimated standard error
#------------------------------------------------------------------------------#
value_OR_se <- function(moOR, data, y, txName, regime) {
			
  #### OLS components
  OLS <- swv_OLS(mo = moOR, data = data, y = y) 

  #### estimator components
  OR <- value_OR_swv(moOR = moOR, 
                     data = data,
                     y = y, 
                     regime = regime, 
                     txName = txName) 

  #### 1,1 Component of Sandwich Estimator

  # OLS contribution
  temp <- OR$dMdB %*% OLS$invdM
  sig11OLS <- temp %*% OLS$MM %*% t(x = temp)

  sig11 <- drop(x = OR$MM + sig11OLS)

  return( c(OR$value, "sigmaHat" = sqrt(x = sig11 / nrow(x = data))) )
  
}
