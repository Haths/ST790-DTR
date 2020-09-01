#------------------------------------------------------------------------------#
# naive estimator of tx effect and its standard error
#
# REQUIRES
#   delta_N()
#
# ASSUMPTIONS
#   tx is binary coded as {0,1}
#
# INPUTS
#  y : a vector containing the outcome of interest
#  A : a vector containing the tx received
#
# RETURNS
#  a list containing
#  deltaHat : the naive estimator for the tx effect
#        EY : the sample mean of the outcome under each tx option
#  sigmaHat : the estimated standard error
#------------------------------------------------------------------------------#
delta_N_se <- function(y, A) {

  # tx effect
  delta <- delta_N(y = y, A = A)

  # estimated tx effect each individual
  delta_i <- y * A / mean(x = A) - y * {1.0 - A} / mean(x = 1.0 - A)

  # variance
  sigmaSq <- mean(x = {delta_i - delta$deltaHat}^2L)

  return( c(delta, "sigmaHat" = sqrt(x = sigmaSq / length(x = y))) ) 

}
