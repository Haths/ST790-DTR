#------------------------------------------------------------------------------#
# tx effect estimator component of the sandwich estimator 
# doubly robust estimator
#
# REQUIRES
#   delta_DR()
#
# ASSUMPTIONS
#   outcome regression model is a linear model
#   propensity score regression model is a logistic model
#   tx is denoted by A and is binary coded as {0,1}
#
# INPUTS:
#  moOR : a modeling object specifying the outcome regression step
#         *** must be a linear model ***
#  moPS : a modeling object specifying the propensity score regression step
#         *** must be a logistic model ***
#  data : a data.frame containing covariates and tx
#         *** tx must be named 'A' and coded 0/1 ***
#     y : a vector containing the outcome of interest
#
# OUTPUTS:
#  a list containing
#  delta : the list returned by delta_DR()
#     MM : M M^T matrix
#   dMdB : matrix of derivatives of M wrt beta
#   dMdG : matrix of derivatives of M wrt gamma
#------------------------------------------------------------------------------#
delta_DR_swv <- function(moOR, moPS, data, y) {

  # estimate treatment effect
  delta <- delta_DR(moOR = moOR, moPS = moPS, data = data, y = y)

  # pi(x;gammaHat)
  p1 <- modelObj::predict(object = delta$fitPS, newdata = data)
  if (is.matrix(x = p1)) p1 <- p1[,ncol(x = p1), drop = TRUE]
          
  # model.matrix for propensity model
  piMM <- stats::model.matrix(object = modelObj::model(object = moPS), 
                              data = data)
          
  A <- data$A
          
  # Q(x,A=0;betaHat)
  data$A <- 0L 
  Q0 <- drop(modelObj::predict(object = delta$fitOR, newdata = data))

  # dQ(x,A=0;betaHat)
  # derivative of a linear model is the model.matrix
  dQ0 <- stats::model.matrix(object = modelObj::model(object = moOR), 
                             data = data)

  # Q(x,A=1;betaHat)
  data$A <- 1L 
  Q1 <- drop(modelObj::predict(object = delta$fitOR, newdata = data))

  # dQ(x,A=1;betaHat)
  # derivative of a linear model is the model.matrix
  dQ1 <- stats::model.matrix(object = modelObj::model(object = moOR), 
                             data = data)

  data$A <- A

  # delta component of M-equation
  mDeltai <- data$A*y/p1 - {1.0 - data$A}*y/{1.0 - p1} - 
             {data$A - p1}*{Q1/p1 + Q0/{1.0 - p1}} - delta$deltaHat

  mmDelta <- mean(x = mDeltai^2)
          
  # derivative of delta component w.r.t. beta
  dMdB <- colMeans(x = {data$A - p1}*{dQ1/p1 + dQ0/{1.0 - p1}})
          
  # derivative of delta component w.r.t. gamma
  dMdG <- - data$A/p1^2*{y - Q1} - {1.0 - data$A}/{1.0 - p1}^2*{y - Q0} 
  dMdG <- colMeans(x = dMdG*p1*{1.0 - p1}*piMM) 
          
  return( list("delta" = delta,
                  "MM" = mmDelta,
                "dMdB" = dMdB,
                "dMdG" = dMdG) )

}
