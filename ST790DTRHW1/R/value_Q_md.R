#------------------------------------------------------------------------------#
# regression based estimator for value of a fixed tx regime with feasible tx 
#   sets at each of the K>1 decision points
#
# ASSUMPTIONS
#   each tx is coded as {0,1}
#
# INPUTS
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
#      fitQ : list of value objects returned by Q-function regressions
#  valueHat : estimated value
#------------------------------------------------------------------------------#
value_Q_md <- function(moQ, data, y, txName, regime) {
  
  # number of decision points
  K <- length(x = moQ)

  #### Q-function Regressions

  qStep <- function(moQ, data, y, txName, regime) {
    # fit Q-function model
    fitObj <- modelObj::fit(object = moQ, data = data, response = y)

    # set tx received to recommended tx
    data[,txName] <- regime

    # predict Q-function for recommended tx
    qk <- modelObj::predict(object = fitObj, newdata = data)

    return (list("fit" = fitObj, "qk" = qk))
  }

  fitQ <- list()
  response <- y

  for (k in K:2L) {
    
    ### consider only individual with more than 1 tx option
    one <- data[,txName[k-1L]] == 1L

    ## kth Q-function regression

    sk <- qStep(moQ = moQ[[ k ]],
                data = data[!one,,drop = FALSE],
                y = response[!one],
                txName = txName[k],
                regime = regime[!one,k])

    fitQ[[ k ]] <- sk$fit
    response[!one] <- sk$qk

  }

  ## k=1 Q-function regression

  s1 <- qStep(moQ = moQ[[ 1L ]],
              data = data,
              y = response,
              txName = txName[1L],
              regime = regime[,1L])

  fitQ[[ 1L ]] <- s1$fit
  response <- s1$qk

  value <- mean(x = response)

  return( list("fitQ" = fitQ, "valueHat" = value) )
}
