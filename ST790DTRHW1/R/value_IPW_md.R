#------------------------------------------------------------------------------#
# IPW estimator for value of a fixed tx regime with binary tx at each of the
#   K>1 decision points
#
# ASSUMPTIONS
#   each tx is binary coded as {0,1}
#
# INPUTS
#  moPS   : list of modeling objects for propensity score regressions
#           kth element corresponds to kth decision point
#  data   : data.frame containing all covariates and txs
#           *** each tx must be coded 0/1 ***
#  y      : outcome of interest
#  txName : vector of tx names (each tx must be coded as 0/1)
#  regime : matrix of 0/1 containing recommended txs
#           kth column corresponds to kth decision point
#
# RETURNS
#  a list containing
#     fitPS : list of value objects returned by propensity score regression
#  valueHat : estimated value
#------------------------------------------------------------------------------#
value_IPW_md <- function(moPS, data, y, txName, regime) {

  # number of decision points
  K <- length(x = moPS)

  # number of patients in data set
  n <- nrow(x = data)

  #### Propensity Score Regressions

  piStep <- function(moPS, data, y, regime) {
    # fit propensity score model
    fitObj <- modelObj::fit(object = moPS, data = data, response = y)

    # predict propensity score
    pk <- modelObj::predict(object = fitObj, newdata = data)
    if (is.matrix(x = pk)) pk <- pk[,ncol(x = pk),drop = TRUE]

    # pull propensity score corresponding to recommended treatment
    pidk <- pk*{regime == 1L} + {1.0 - pk}*{regime == 0L}

    return (list("fit" = fitObj, "pidk" = pidk))
  }

  fitPS <- list()
  pid <- matrix(data = 1.0, nrow = n, ncol = K)

  for (k in K:2L) {

    # use only individuals with more than 1 tx option
    one <- data[,txName[k-1L]] == 1L

    ## kth propensity score regression

    sk <- piStep(moPS = moPS[[ k ]], 
                 data = data[!one,,drop = FALSE],  
                 y = data[!one,txName[k]],
                 regime = regime[!one,k])

    fitPS[[ k ]] <- sk$fit
    pid[!one,k] <- sk$pidk

  }

  ## k=1 propensity score regression

  s1 <- piStep(moPS = moPS[[ 1L ]], 
               data = data,  
               y = data[,txName[1L]], 
               regime = regime[,1L])

  fitPS[[ 1L ]] <- s1$fit
  pid[,1L] <- s1$pidk

  # indicator of all txs received being those recommended
  cee <- apply(X = {data[,txName] == regime}, MARGIN = 1L, FUN = prod)

  # P(C == 1|H)
  pid <- apply(X = pid, MARGIN = 1L, FUN = prod)

  #### Value Estimate

  value <- mean(x = cee * y / pid)

  return( list("fitPS" = fitPS, "valueHat" = value) )
}
