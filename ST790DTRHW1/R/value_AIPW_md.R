#------------------------------------------------------------------------------#
# AIPW estimator for value of a fixed tx regime with binary tx at each of the
#   K>1 decision points
#
# ASSUMPTIONS
#   each tx is binary coded as {0,1}
#
# INPUTS
#  moPS   : list of modeling objects for propensity score regressions
#           must be logistic models; kth element corresponds to kth decision
#           point
#  moQ    : list of modeling objects for Q-function regressions
#           kth element corresponds to kth decision point
#  data   : data.frame containing all covariates and txs
#  y      : outcome of interest
#  txName : vector of tx names (each tx must be coded as 0/1)
#  regime : matrix of 0/1 containing recommended txs
#           kth column corresponds to kth decision point
#
# RETURNS
#  a list containing
#     fitPS : list of value objects returned by propensity score regression
#      fitQ : list of value objects returned by Q-function regressions
#  valueHat : estimated value
#------------------------------------------------------------------------------#
value_AIPW_md <- function(moPS, moQ, data, y, txName, regime) {

  # number of decision points
  K <- length(x = moPS)

  # number of patients in data set
  n <- nrow(x = data)

  #### Propensity Score and Q-function Regressions

  qStep <- function(moQ, data, y, txName, regime) {
    # fit Q-function model
    fitObj <- modelObj::fit(object = moQ, data = data, response = y)

    # set tx received to recommended tx
    data[,txName] <- regime

    # predict Q-function for recommended tx
    qk <- modelObj::predict(object = fitObj, newdata = data)

    return (list("fit" = fitObj, "qk" = qk))
  }

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

  fitQ <- list()
  res <- matrix(data = 0.0, nrow = n, ncol = K+1L)

  fitPS <- list()
  pid <- matrix(data = 1.0, nrow = n, ncol = K)


  res[,K+1L] <- y

  for (k in K:2L) {

    # use only individuals with more than 1 tx option
    one <- data[,txName[k-1L]] == 1L

    ## kth Q-function regression

    sk <- qStep(moQ = moQ[[ k ]],
                data = data[!one,,drop = FALSE],
                y = res[!one,k+1L],
                txName = txName[k],
                regime = regime[!one,k])

    fitQ[[ k ]] <- sk$fit

    res[,k] <- res[,k+1L]
    res[!one,k] <- sk$qk

    ## kth propensity score regression

    sk <- piStep(moPS = moPS[[ k ]], 
                 data = data[!one,,drop = FALSE],  
                 y = data[!one,txName[k]],
                 regime = regime[!one,k])

    fitPS[[ k ]] <- sk$fit
    pid[!one,k] <- sk$pidk

  }

  # k=1 Q-function regression
  s1 <- qStep(moQ = moQ[[ 1L ]],
              data = data,
              y = res[,2L],
              txName = txName[1L],
              regime = regime[,1L])

  fitQ[[ 1L ]] <- s1$fit
  res[,1L] <- s1$qk

  # k=1 propensity score regression
  s1 <- piStep(moPS = moPS[[ 1L ]], 
               data = data,  
               y = data[,txName[1L]], 
               regime = regime[,1L])

  fitPS[[ 1L ]] <- s1$fit
  pid[,1L] <- s1$pidk

  # indicators of 1:k txs received being those recommended
  cee <- t(x = apply(X = {data[,txName] == regime}, 
                     MARGIN = 1L,  
                     FUN = cumprod))

  # P(Ck == 1|Hk)
  pid <- t(x = apply(X = pid, MARGIN = 1L, FUN = cumprod))

  #### Value Estimate

  value <- cee[,K] * y / pid[,K]
  value <- value + {1.0 - cee[,1L] / pid[,1L]} * res[,1L]

  for (k in 2L:K) {
    value <- value + {cee[,k-1L] / pid[,k-1L] - cee[,k] / pid[,k]}*res[,k]
  }

  value <- mean(x = value)

  return( list("fitPS" = fitPS, "fitQ" = fitQ, "valueHat" = value) )
}
