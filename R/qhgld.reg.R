#' @import GLDEX
#' @export
#' @title Predict quantiles of a Hurdle Generalized Lambda regression
#'
#' @description Quantile function of a Hurdle Generalized Lambda Regression.
#'
#' @details Given an object of class \link[HGLD]{reg.hgld}, a percentile and new values for the covariates,
#'  it returns the quantile of the fitted HGLD, given the covariates. The fitted HGLD must be non-negative. The
#'  contrast on the \emph{newvalues} data frame must be the same of the \emph{data} used to fit the \link[HGLD]{reg.hgld} regression.
#'
#' @param p Vector of probabilities.
#' @param reg An object of class \link[HGLD]{reg.hgld}.
#' @param newvalues A data frame with the new values of the covariates. Column names must match the ones given in formulas \emph{loc.formula} and \emph{zero.formula} of
#' the \link[HGLD]{reg.hgld} object.
#' @param l0 Whether to return the lambda0 of each profile.
#' @param location Whether to return the location of each profile.
#' @return The quantile of the new values, according to the fitted HGLD regression.
#' @examples
#' set.seed(100)
#' tmp <- na.omit(healthcare)
#' data <- tmp[sample(1:nrow(tmp),50),]
#' formula <- log_expense ~ age + sex + log_previous_expense
#' reg <- suppressWarnings(reg.hgld(data = data,zero.formula = formula,loc.formula = formula,
#'                         full = FALSE,param = "rs"))
#' newvalues <- tmp[sample(1:nrow(tmp),5),c(2,3,8)]
#' qhgld.reg(p = seq(0.05,0.95,0.1),reg,newvalues)
#'
#' @references Marcondes, D.; Peixoto, C.; Maia, A. C.; A Survey of a Hurdle Model for Heavy-Tailed Data Based on the Generalized Lambda Distribution. (2017) \emph{arxiv1712.02183}

qhgld.reg <- function(p,reg,newvalues,l0 = TRUE,location = TRUE){
  # Estimated coeffcients
  if(reg$full){
    lambda <- reg$GLDreg[[1]]$`Estimated parameters`
    lambda <- lambda[(length(lambda)-3):length(lambda)]
    beta <- reg$GLDreg[[1]]$`Estimated parameters`
    beta <- beta[1:(length(beta)-4)]
  }
  else{
    lambda <- reg$GLDreg$`Estimated parameters`
    lambda <- lambda[(length(lambda)-3):length(lambda)]
    beta <- reg$GLDreg$`Estimated parameters`
    beta <- beta[1:(length(beta)-4)]
  }
  gama <- reg$gamlss$mu.coefficients
  param <- reg$param

  # Predictors
  newvalues[[all.vars(reg$loc.formula)[1]]] <- rep(1,nrow(newvalues))
  W <- model.matrix(lm(reg$loc.formula,data = newvalues))
  Z <- model.matrix(lm(reg$zero.formula,data = newvalues))
  local <- W %*% beta
  odds <- exp(Z %*% gama)
  lambda0 <- odds/(1+odds)
  p <- data.frame("p" = p,"order" = rank(p))
  p <- p[order(p$order),]

  # Quantiles
  quantiles <- matrix(nrow = length(lambda0),ncol = nrow(p))
  for(prof in 1:length(lambda0)){
    pzero <- pgl(q = -local[prof],lambda1 = lambda,param = param)
    if(pzero > 0){
      cat(paste("The distribution fitted to profile",prof,"is not positive!\n"))
    }
    else{
      quantiles[prof,p$p <= lambda0[prof]] <- 0
      quantiles[prof,p$p > lambda0[prof]] <- qgl(p = (p$p[p$p > lambda0[prof]] - lambda0[prof])/(1 - lambda0[prof]),lambda1 = lambda,param = param) + local[prof]
    }
  }
  colnames(quantiles) <- p$p
  rownames(quantiles) <- rownames(newvalues)
  quantiles <- quantiles[,p$order]
  if(location){
    n <- colnames(quantiles)
    quantiles <- data.frame("Location" = local,quantiles)
    colnames(quantiles) <- c("Location",n)
  }
  if(l0){
    n <- colnames(quantiles)
    quantiles <- data.frame("lambda0" = lambda0,quantiles)
    colnames(quantiles) <- c("lambda0",n)
  }
  return(quantiles)
}
