#------------------------------------------------------------------------------#
# component of sandwich estimator for IPW estimator of value of a fixed regime
#
# REQUIRES
#   value_IPW()
#
# ASSUMPTIONS
#   propensity score model is a logistic model
#   tx is binary coded as {0,1}
#
# INPUTS
#  moPS   : modeling object for propensity score regression
#           *** must be a logistic model ***
#  data   : data.frame containing baseline covariates and tx
#           *** tx must be coded 0/1 ***
#  y      : outcome of interest
#  txName : tx column in data (tx must be coded as 0/1)
#  regime : 0/1 vector containing recommended tx
#
# RETURNS
#  a list containing
#  value : list returned by value_IPW()
#     MM : M M^T matrix
#   dMdG : matrix of derivatives of M wrt gamma
#------------------------------------------------------------------------------#
value_IPW_swv <- function(moPS, data, y, txName, regime) {

  # estimate value
  value <- value_IPW(moPS = moPS, 
                     data = data,  
                     y = y,  
                     regime = regime,  
                     txName = txName)
  
  # pi(x; gamma)
  p1 <- modelObj::predict(object = value$fitPS, newdata = data) 
  if (is.matrix(x = p1)) p1 <- p1[,ncol(x = p1), drop = TRUE]
  
  # model.matrix for propensity model
  piMM <- stats::model.matrix(object = modelObj::model(object = moPS), 
                              data = data)
  
  # propensity for receiving recommended tx
  pid <- p1*{regime == 1L} + {1.0-p1}*{regime == 0L}
  
  # indicator of tx received = regime
  Cd <- regime == data[,txName]
  
  # value component of M MT
  mmValue <- mean(x = {Cd * y / pid - value$valueHat}^2)
  
  # derivative w.r.t. gamma
  dMdG <- colMeans(x = Cd * y / pid * (-regime + p1) * piMM)
  
  return( list("value" = value, "MM" = mmValue, "dMdG" = dMdG) )

}
