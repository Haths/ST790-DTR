#------------------------------------------------------------------------------#
# Alternative IPW estimator for value of a fixed binary tx regime
#
# ASSUMPTIONS
#   tx is binary coded as {0,1}
#
# INPUTS
#    moPS : modeling object specifying propensity score regression
#    data : data.frame containing baseline covariates and tx
#           *** tx must be coded 0/1 ***
#       y : outcome of interest
#  txName : tx column in data (tx must be coded as 0/1)
#  regime : 0/1 vector containing recommended tx
#
# RETURNS
#  a list containing
#     fitPS : value object returned by propensity score regression
#        EY : sample average outcome for each tx
#  valueHat : estimated value
#------------------------------------------------------------------------------#
value_A_IPW <- function(moPS, data, y, txName, regime) {
  
  #### Propensity Score Regression
  fitPS <- modelObj::fit(object = moPS, data = data, response = data[,txName])
  
  # estimate propensity score
  p1 <- modelObj::predict(object = fitPS, newdata = data)
  if (is.matrix(x = p1)) p1 <- p1[,ncol(x = p1), drop = TRUE]
  
  #### Value Estimate
  
  EY <- array(data = 0.0, dim = 2L, dimnames = list(c("0","1")))
  sub1 <- {regime == data[,txName]} & {data[,txName] == 1L}
  sub0 <- {regime == data[,txName]} & {data[,txName] == 0L}
  EY[2L] <- sum(x = sub1 * y / p1)
  EY[1L] <- sum(x = sub0 * y / {1.0 - p1})
  Weight[2L] <- sum(x = sub1  / p1)
  Weight[1L] <- sum(x = sub0  / {1.0 - p1})
  
  value <- sum(EY)/sum(Weight)
  
  return( list("fitPS" = fitPS, "EY" = EY, "valueHat" = value) )
}