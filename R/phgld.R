#' @import GLDEX
#' @export
#' @title The Hurdle Generalized Lambda Distribution Family
#'
#' @description Probability function of the Hurdle Generalized Lambda Distribution.
#'
#' @details If the parametrization of the RS or FMKL GLD is not acceptable, the function returns NA. This function is based on the \link[GLDEX]{GLDEX} package.
#'
#' @param q Vector of quantiles.
#' @param mixture Whether to give the density of a mixture HGLD.
#' @param lambda1 A vector of length 5 with the five parameters of the HGLD or the first HGLD if mixture = TRUE.
#' @param lambda2 A vector of length 4 with the four parameters of the second HGLD.
#' @param prob The cluster parameter for the mixture HGLD.
#' @param param "fmkl" or "rs".
#' @param inverse.eps Accuracy of calculation for the numerical determination of F(x), defaults to 1e-8.
#' @param max.iterations Maximum number of iterations in the numerical determination of F(x), defaults to 500.
#' @return The probability function of the HGLD.
#' @examples
#' library(ggplot2)
#' set.seed(100)
#' data <- healthcare[sample(1:nrow(healthcare),100),]
#' fit <- fit.hgld(data$log_expense)
#' {ggplot(data.frame(x = data$log_expense),aes(x = x)) +
#' stat_function(fun = function(x) phgld(q = x,lambda1 = fit$par[,2],param = "rs"))}
#'
#' #mixture
#' set.seed(100)
#' data <- c(rcauchy(50,location = 10),rep(0,30),rcauchy(50))
#' fit <- fit.hgld(data,TRUE)
#' {ggplot(data.frame(x = data[data!=0]),aes(x = x)) +
#' stat_function(fun = function(x) phgld(q = x,mixture = TRUE,lambda1 = fit$par[1:5,2],
#' lambda2 = fit$par[6:9,2],prob = fit$par[10,2],param = "rs"))}
#'
#' @references Marcondes, D.; Peixoto, C.; Maia, A. C.; A Survey of a Hurdle Model for Heavy-Tailed Data Based on the Generalized Lambda Distribution. (2017) \emph{arxiv1712.02183}
#' @references Su, S.; Fitting Single and Mixture of Generalized Lambda Distributions to Data via Discretized and Maximum Likelihood Methods: GLDEX in R. (2007), Journal of Statistical Software: *21* 9.

phgld <- function(q,mixture = FALSE,lambda1,lambda2 = NULL,prob = NULL,param = "fmkl",inverse.eps = 1e-08,max.iterations = 500){
  if(mixture == TRUE){
    r <- (1 - lambda1[1]) *
      (prob * pgl(q = q,lambda1 = lambda1[-1],param = param,inverse.eps = inverse.eps,max.iterations = max.iterations) +
         (1 - prob) * pgl(q = q,lambda1 = lambda2,param = param,inverse.eps = inverse.eps,max.iterations = max.iterations))
    r[q >= 0] <- r[q >= 0] + lambda1[1]
    return(r)
  }
  else{
    r <- (1 - lambda1[1]) * pgl(q = q,lambda1 = lambda1[-1],param = param,inverse.eps = inverse.eps,max.iterations = max.iterations)
    r[q >= 0] <- r[q >= 0] + lambda1[1]
    return(r)
  }
}
