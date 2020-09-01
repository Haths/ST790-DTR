#------------------------------------------------------------------------------#
# naive estimator of tx effect
#
# ASSUMPTIONS
#   tx is binary coded as 0,1
#
# INPUTS
#  y : a vector containing the outcome of interest
#  A : a vector containing the tx received
#
# RETURNS
#  a list containing
#        EY : the sample mean of the outcome under each tx option
#  deltaHat : the naive estimator for the tx effect
#------------------------------------------------------------------------------#
delta_N <- function(y, A) {

  #### Average Responses

  # aggregate data for mean in each tx
  EY <- stats::aggregate(x = y, by = list(A), FUN = mean)

  # convert to named vector
  EY <- array(data = EY[,2L], dimnames = list(EY[,1L]))

  #### Treatment Effect

  delta <- unname(obj = EY[2L] - EY[1L])

  return( list("EY" = EY, "deltaHat" = delta) )
}
