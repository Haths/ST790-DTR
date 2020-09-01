#------------------------------------------------------------------------------#
# tx effect estimator component of the sandwich estimator IPW estimator
#
# REQUIRES
#   delta_IPW()
#
# ASSUMPTIONS
#   propensity score regression model is a logistic model
#   tx is denoted by A and is binary coded as {0,1}
#
# INPUTS:
#  moPS : a modeling object specifying the propensity score regression step
#         *** must be a logistic model ***
#  data : a data.frame containing covariates and tx
#         *** tx must be named 'A' and coded 0/1 ***
#     y : a vector containing the outcome of interest
#
# OUTPUTS:
#  a list containing
#  delta : the list returned by delta_IPW()
#     MM : M M^T matrix
#   dMdG : matrix of derivatives of M wrt gamma
#------------------------------------------------------------------------------#
delta_IPW_swv <- function(moPS, data, y) {

  # estimate treatment effect
  delta <- delta_IPW(moPS = moPS, data = data, y = y)

  # pi(x; gamma)
  p1 <- modelObj::predict(object = delta$fitPS, newdata = data)
  if (is.matrix(x = p1)) p1 <- p1[,ncol(x = p1), drop = TRUE]

  # model.matrix for propensity model
  piMM <- stats::model.matrix(object = modelObj::model(object = moPS), 
                              data = data)

  # delta component of M MT
  mmDelta <- mean(x = {data$A * y / p1 - {1.0-data$A} * y / {1.0-p1} - 
                       delta$deltaHat}^2)

  # derivative w.r.t. gamma
  dMdG <- -{ data$A * y * {1.0-p1} / p1 +
            {1.0-data$A} * y * p1 / {1.0-p1} } * piMM
  dMdG <- colMeans(x = dMdG)

  return( list("delta" = delta,
                  "MM" = mmDelta,
                "dMdG" = dMdG) )

}
