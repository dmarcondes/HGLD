#' @import GLDEX
#' @export
#' @title The Hurdle Generalized Lambda Distribution Family
#'
#' @description Accumulated probability function of the Hurdle Generalized Lambda Distribution.
#'
#' @details If the parametrization of the RS or FMKL GLD is not acceptable, the function returns NA. This function is based on
#' the \link[GLDEX]{GLDEX} package.
#'
#' @param q Vector of quantiles.
#' @param mixture Whether to give the accumulated probability function of a mixture of HGLDs.
#' @param lambda1 A vector of length 5 with the five parameters of the HGLD, or of the first HGLD if mixture = TRUE.
#' @param lambda2 A vector of length 4 with the four parameters of the second HGLD if mixture = TRUE.
#' @param prob The cluster parameter for the mixture of HGLDs.
#' @param param "fmkl" or "rs".
#' @param inverse.eps Accuracy of calculation for the numerical determination of F(x), defaults to 1e-8.
#' @param max.iterations Maximum number of iterations in the numerical determination of F(x), defaults to 500.
#' @return The accumulated probability function of the HGLD.
#' @examples
#' library(ggplot2)
#' {ggplot(data.frame(x = seq(2,8,0.25)),aes(x = x)) +
#' stat_function(fun = function(x) phgld(q = x,lambda1 = c(0.54,4.561,0.0191,0.00925,0.0217),
#'                                       param = "rs"))}
#'
#' #mixture
#' lambda1 <- c(0.230,0.3514,-0.4472,-0.3748,-0.3108)
#' lambda2 <- c(9.624,-1.227,-0.6290,-0.8515)
#' {ggplot(data.frame(x = seq(-5,15,0.25)),aes(x = x)) +
#' stat_function(fun = function(x) phgld(q = x,mixture = TRUE,lambda1 = lambda1,
#' lambda2 = lambda2,prob = 0.479,param = "rs"))}
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
