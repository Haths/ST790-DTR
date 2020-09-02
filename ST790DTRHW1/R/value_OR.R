#------------------------------------------------------------------------------#
# outcome regression estimator for value of a fixed binary tx regime
#
# ASSUMPTIONS
#   tx is binary coded as {0,1}
#
# INPUTS
#    moOR : modeling object specifying outcome regression
#    data : data.frame containing baseline covariates and tx received
#           *** tx must be coded 0/1 ***
#       y : outcome of interest
#  txName : tx column in data (tx must be coded as 0/1)
#  regime : 0/1 vector containing recommended tx
#
# RETURNS
#  a list containing
#     fitOR : value object returned by outcome regression
#        EY : contribution to value from each tx
#  valueHat : estimated value
#------------------------------------------------------------------------------#
value_OR <- function(moOR, data, y, txName, regime) {
  
  #### Outcome Regression 

  fitOR <- modelObj::fit(object = moOR, data = data, response = y)

  # estimated Q-function when all A=d
  data[,txName] <- regime 
  Qd <- drop(x = modelObj::predict(object = fitOR, newdata = data))

  #### Value of regime

  EY <- array(data = 0.0, dim = 2L, dimnames = list(c("0","1")))
  EY[2L] <- mean(x = Qd*{regime == 1L})
  EY[1L] <- mean(x = Qd*{regime == 0L})

  value <- sum(EY)

  return( list("fitOR" = fitOR, "EY" = EY, "valueHat" = value) )
}
