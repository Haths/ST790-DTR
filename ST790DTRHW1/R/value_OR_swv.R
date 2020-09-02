#------------------------------------------------------------------------------#
# component of sandwich estimator for outcome regression estimator of value of
# a fixed regime
#
# REQUIRES
#   value_OR()
#
# ASSUMPTIONS
#   outcome regression model is a linear model
#   tx is binary coded as {0,1}
#
# INPUTS
#  moOR   : modeling object for outcome regression
#           *** must be a linear model ***
#  data   : data.frame containing baseline covariates and tx
#           *** tx must be coded 0/1 ***
#  y      : outcome of interest
#  txName : tx column in data (tx must be coded as 0/1)
#  regime : 0/1 vector containing recommended tx
#
# RETURNS
#  a list containing
#  value : list returned by value_OR()
#     MM : M M^T matrix
#   dMdB : matrix of derivatives of M wrt beta
#------------------------------------------------------------------------------#
value_OR_swv <- function(moOR, data, y, txName, regime) {

  # estimate value
  value <- value_OR(moOR = moOR, 
                    data = data,  
                    y = y,  
                    regime = regime,  
                    txName = txName)
			
  # Q(H,d;betaHat)
  data[,txName] <- regime 
  Qd <- drop(x = modelObj::predict(object = value$fit, newdata = data))
  
  # derivative of Q(H,0;betaHat)
  dQd <- stats::model.matrix(object = modelObj::model(object = moOR), 
                             data = data)
  
  # value component of M MT
  mmt <- mean(x = {Qd - value$valueHat}^2)

  # derivative of value component w.r.t. beta
  dMdB <- colMeans(x = dQd) 
  
  return( list("value" = value, "MM" = mmt, "dMdB" = dMdB) )

}
