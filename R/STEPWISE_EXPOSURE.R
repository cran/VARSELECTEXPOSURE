#' Performs deviance-based stepwise variable selection in logistic regression with an exposure
#'
#' Returns the estimated Average Treatment Effect and estimated Relative Treatment Effect calculated by the optimal model chosen via stepwise selection including an exposure variable.
#' @param Data Data frame containing outcome variable (Y), exposure variable (E), and candidate covariates.
#' @importFrom stats glm pchisq
#' @importFrom utils head
#' @return List containing (1) the estimated Average Treatment Effect, (2) estimated Relative Treatment Effect, (3) summary of the selected model, and (4) the first 6 rows of the data frame containing stepwise-selected covariates.
#' @examples
#' ###Generate data with n rows and p covariates, can be any number but we'll choose 1000 rows
#' ###and 10 covariates for this example
#' set.seed(3)
#'
#' p = 7
#' n = 750
#' beta0 = rnorm(1, mean = 0, sd = 1)
#' betaE = rnorm(1, mean = 0, sd = 1)
#' beta0_E = rnorm(1, mean = 0, sd = 1)
#' betaX_E = c()
#' betaX_Y = c()
#' Y = rep(NA, n)
#' E = rep(NA, n)
#' pi0 = rep(NA, n)
#' pi1 = rep(NA, n)
#' data = data.frame(cbind(Y, E, pi0, pi1))
#' j = round(runif(1, 0, p))
#' for(i in 1:p){
#'   betaX_Y[i] = rnorm(1, mean = 0, sd = 0.5)
#'   betaX_E[i] = rnorm(1, mean = 0, sd = 0.5)
#' }
#' zeros = sample(1:p, j, replace =  FALSE)
#' betaX_Y[zeros] = 0
#' betaX_E[zeros] = 0
#' mu = 0
#' sigma = 1
#' for(i in 1:p){
#'   covar = rnorm(n, 0, 1)
#'   data[,i+4] = covar
#'   names(data)[i+4] = paste("X", i, sep = "")
#' }
#' for(i in 1:n){
#'   p.event_E = beta0_E + sum(betaX_E*data[i,5:(p+4)])
#'   pi1_E = exp(p.event_E)/(1+exp(p.event_E))
#'   data[i,2] = rbinom(1, 1, prob = pi1_E)
#' }
#' for(i in 1:n){
#'   p.event = beta0 + betaE + sum(betaX_Y*data[i,5:(p+4)])
#'   p.noevent = beta0 + sum(betaX_Y*data[i,5:(p+4)])
#'   pi0 = exp(p.noevent)/(1+exp(p.noevent))
#'   pi1 = exp(p.event)/(1+exp(p.event))
#'   data[i,3] = pi0
#'   data[i,4] = pi1
#'   if(data[i,2] == 1){
#'     data[i,1] = rbinom(1, 1, prob = pi1)
#'   }else{
#'     data[i,1] = rbinom(1, 1, prob = pi0)
#'   }
#' }
#' for(i in 1:n){
#'   p.event_E = beta0_E + sum(betaX_E*data[i,5:(p+4)])
#'   pi1_E = exp(p.event_E)/(1+exp(p.event_E))
#'   data[i,2] = rbinom(1, 1, prob = pi1_E)
#' }
#'
#' ###Raw data includes pi0 and pi1 columns used to fill Y and E, so to test
#' ###the function we'll remove these
#'
#' testdata = data[,-c(3,4)]
#'
#' STEPWISE_EXPOSURE(testdata)
#' @export
STEPWISE_EXPOSURE = function(Data){
  p = ncol(Data) - 2
  n = nrow(Data)
  count = 0
  VAR.NUM = 3
  VAR.INDEX = seq(1, p, by = 1)
  DEV.MOD = glm(Y~1+E, data = Data, family = "binomial")
  D0 = DEV.MOD$deviance
  VAR.INDEX = sort(VAR.INDEX)
  STEPWISE.DATA = Data[,c(1,2)]
  while(count != 2){
    count = 0
    DEVS = c()
    for(VAR in VAR.INDEX){
      test.data = cbind(STEPWISE.DATA, Data[,VAR+2])
      test.mod = glm(Y~., data = test.data, family = "binomial")
      DEVS = append(DEVS, test.mod$deviance)
    }
    new.index = VAR.INDEX[which.min(DEVS)]
    D1 = min(DEVS)
    TEST = abs(D1-D0)
    pval = 1-pchisq(TEST, 1)
    if(pval < 0.05){
      STEPWISE.DATA = cbind(STEPWISE.DATA, Data[,(new.index+2)])
      colnames(STEPWISE.DATA)[VAR.NUM] = names(Data)[VAR.NUM]
      VAR.NUM = VAR.NUM + 1
      D0 = D1
      VAR.INDEX = VAR.INDEX[-which.min(DEVS)]
      count = 0
    }else{
      count = count + 1
    }


    DEVS = c()
    for(VAR in 1:(VAR.NUM-3)){
      test.data = STEPWISE.DATA[,-(VAR+2)]
      test.mod = glm(Y~., data = test.data, family = "binomial")
      DEVS = append(DEVS, test.mod$deviance)
    }
    new.index = which.min(DEVS)
    D1 = min(DEVS)
    TEST = abs(D1-D0)
    pval = 1-pchisq(TEST, 1)
    if(pval > 0.05){
      STEPWISE.DATA = STEPWISE.DATA[,-(new.index+2)]
      D0 = D1
      VAR.NUM = VAR.NUM - 1
      VAR.INDEX = append(VAR.INDEX, new.index)
      count = 0
    }else{
      count = count + 1
    }
  }
  MOD.LENGTH = length(colnames(STEPWISE.DATA))
  LOGIT.MOD = glm(Y~., data = STEPWISE.DATA, family = "binomial")
  BETA.EST = LOGIT.MOD$coefficients
  pi0.E = c()
  pi1.E = c()
  for(j in 1:n){
    p.event = as.numeric(BETA.EST[1] + BETA.EST[2] + sum(BETA.EST[3:MOD.LENGTH]*STEPWISE.DATA[j,3:MOD.LENGTH]))
    p.noevent = as.numeric(BETA.EST[1] + sum(BETA.EST[3:MOD.LENGTH]*STEPWISE.DATA[j,3:MOD.LENGTH]))
    pi0 = exp(p.noevent)/(1+exp(p.noevent))
    pi1 = exp(p.event)/(1+exp(p.event))
    pi0.E[j] = pi0
    pi1.E[j] = pi1
  }
  RESULT = vector(mode = "list")
  RESULT$ATE = sum(pi1.E - pi0.E)/n
  RESULT$RTE = sum(pi1.E/pi0.E)/n
  RESULT$DATA = head(STEPWISE.DATA)
  RESULT$MOD = summary(LOGIT.MOD)
  return(RESULT)
}
