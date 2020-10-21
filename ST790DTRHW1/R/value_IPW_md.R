#------------------------------------------------------------------------------#
# IPW estimator for value of a fixed tx regime with binary tx at each of the
#   K>1 decision points
#
# ASSUMPTIONS
#   each of the tx is denoted by Ak k=(1,...,K) and is binary coded as {0,1}
#
# INPUTS
#  moPS   : list of modeling objects for propensity score regressions
#           must be logistic models; kth element corresponds to kth decision
#           point
#  data   : data.frame containing all covariates and txs
#           txs must be named 'Ak' (k=1,...,K) and coded 0/1
#  y      : outcome of interest
#  regime : matrix of 0/1 containing recommended txs
#           kth column corresponds to kth decision point
#
# RETURNS
#  a list containing
#     fitPS : list of value objects returned by propensity score regression
#  valueHat : estimated value
#------------------------------------------------------------------------------#
value_IPW_md <- function(moPS, data, y, regime) {

  # number of decision points
  K <- length(x = moPS)

  # number of patients in data set
  n <- nrow(x = data)

  # tx names in data set
  nms <- paste0('A',1L:K)

  #### Propensity Score Regression

  fitPS <- list()
  pid <- rep(x = 1.0, times = n)

  for (k in 1L:K) {

    # fit propensity score
    fitPS[[ k ]] <- modelObj::fit(object = moPS[[ k ]], 
                                  data = data, 
                                  response = data[,nms[k]])

    # product of estimated propensity scores
    pk <- modelObj::predict(object = fitPS[[ k ]], newdata = data)
    if (is.matrix(x = pk)) pk <- pk[,ncol(x = pk), drop = TRUE]

    pid <- pid * {pk*{regime[,k] == 1L} + {1.0 - pk}*{regime[,k] == 0L}}

  }

  #### Value Estimate

  # indicator of tx received was that dictated by the regime
  Cd <- rowSums(x = regime == data[,nms]) == K

  value <- mean(x = Cd * y / pid)

  return( list(   "fitPS" = fitPS,
               "valueHat" = value) )
}
