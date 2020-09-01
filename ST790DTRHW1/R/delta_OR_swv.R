#------------------------------------------------------------------------------#
# tx effect estimator component of the sandwich estimator 
# outcome regression estimator
#
# REQUIRES
#   delta_OR()
#
# ASSUMPTIONS
#   outcome regression model is a linear model
#   tx is denoted by A and is binary coded as {0,1}
#
# INPUTS:
#  moOR : a modeling object specifying the outcome regression step
#         *** must be a linear model ***
#  data : a data.frame containing covariates and tx
#         *** tx must be named 'A' and coded 0/1 ***
#     y : a vector containing the outcome of interest
#
# OUTPUTS:
#  a list containing
#  delta : the list returned by delta_OR()
#     MM : M M^T matrix
#   dMdB : matrix of derivatives of M wrt beta
#------------------------------------------------------------------------------#
delta_OR_swv <- function(moOR, data, y) {

  # estimate tx effect
  delta <- delta_OR(moOR = moOR, data = data, y = y)

  # Q(X,0;betaHat)
  data$A <- 0L 
  Q0 <- drop(modelObj::predict(object = delta$fitOR, newdata = data))

  # derivative of Q(X,0;betaHat)
  dQ0 <- stats::model.matrix(object = modelObj::model(object = moOR), 
                             data = data)

  # Q(X,1;betaHat)
  data$A <- 1L 
  Q1 <- drop(modelObj::predict(object = delta$fitOR, newdata = data))

  # derivative of Q(X,1;betaHat)
  dQ1 <- stats::model.matrix(object = modelObj::model(object = moOR), 
                             data = data)

  # delta component of M MT
  mmDelta <- mean(x = {Q1 - Q0 - delta$deltaHat}^2)

  # derivative of delta component w.r.t beta
  dMdB <- colMeans(x = dQ1 -dQ0) 

  return( list("delta" = delta, "MM" = mmDelta, "dMdB" = dMdB) )

}
