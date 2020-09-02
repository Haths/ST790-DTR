#------------------------------------------------------------------------------#
# value and standard error for the IPW estimator of the value of a fixed regime
# using the sandwich estimator of variance
#
# REQUIRES
#   swv_ML() and value_IPW_swv()
#
# ASSUMPTIONS
#   tx is binary coded as {0,1}
#   moPS is a logistic model
#
# INPUTS
#    moPS : modeling object specifying propensity score regression
#           *** must be a logistic model ***
#    data : data.frame containing covariates and tx
#           *** tx must be coded 0/1 ***
#       y : vector of outcome of interest
#  txName : tx column in data (tx must be coded as 0/1)
#  regime : 0/1 vector containing recommended tx
#
# RETURNS
#  a list containing
#     fitPS : modelObjFit object for propensity score regression
#        EY : sample average outcome for each tx
#  valueHat : estimated value
#  sigmaHat : estimated standard error
#------------------------------------------------------------------------------#
value_IPW_se <- function(moPS, data, y, txName, regime) {

  # obtain ML components
  ML <- swv_ML(mo = moPS, data = data, y = data[,txName]) 

  # obtain IPW value components
  IPW <- value_IPW_swv(moPS = moPS, 
                       data = data,  
                       y = y,  
                       regime = regime,  
                       txName = txName) 

  # calculate 1,1 component
  temp <- IPW$dMdG %*% ML$invdM
  sig11 <- drop(x = IPW$MM + temp %*% ML$MM %*% t(x = temp))

  return( c(IPW$value, "sigmaHat" = sqrt(x = sig11/ nrow(x = data))) )
}
