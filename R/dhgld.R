#' @import GLDEX
#' @export
#' @title The Hurdle Generalized Lambda Distribution Family
#'
#' @description Density of the Hurdle Generalized Lambda Distribution.
#'
#' @details If the parametrization of the RS or FMKL HGLD is not acceptable, the function returns NA. This function is based on the
#' \link[GLDEX]{GLDEX} package.
#'
#' @param x Vector of data.
#' @param mixture Whether to give the density of a mixture of HGLDs.
#' @param lambda1 A vector of length 5 with the five parameters of the HGLD, or of the first HGLD if mixture = TRUE.
#' @param lambda2 A vector of length 4 with the four parameters of the second HGLD if mixture = TRUE.
#' @param prob The cluster parameter (probability) for the mixture of HGLDs.
#' @param param "fmkl" or "rs".
#' @param inverse.eps Accuracy of calculation for the numerical determination of F(x), defaults to 1e-8.
#' @param max.iterations Maximum number of iterations in the numerical determination of F(x), defaults to 500.
#' @return The probability density of the continuous part of the HGLD.
#' @examples
#' library(ggplot2)
#' {ggplot(data.frame(x = seq(-3,1,0.01)),aes(x = x)) +
#' stat_function(fun = function(x) dhgld(x = x,
#'                                      lambda1 = c(0.48,0.509,-0.000369,-0.0002483,-0,0003916),
#'                                      param = "rs"))}
#'
#' #mixture
#' lambda1 <- c(0.230,0.3514,-0.4472,-0.374,-0.3108)
#' lambda2 <- c(9.624,-1.227,-0.629,-0.8515)
#' {ggplot(data.frame(x = seq(-10,20,0.25)),aes(x = x)) +
#' stat_function(fun = function(x) dhgld(x = x,mixture = TRUE,lambda1 = lambda1,
#' lambda2 = lambda2,prob = 0.47954,param = "rs"))}
#'
#' @references Marcondes, D.; Peixoto, C.; Maia, A. C.; A Survey of a Hurdle Model for Heavy-Tailed Data Based on the Generalized Lambda Distribution. (2017) \emph{arxiv1712.02183}
#' @references Su, S.; Fitting Single and Mixture of Generalized Lambda Distributions to Data via Discretized and Maximum Likelihood Methods: GLDEX in R. (2007), Journal of Statistical Software: *21* 9.

dhgld <- function(x,mixture = FALSE,lambda1,lambda2 = NULL,prob = NULL,param = "fmkl",inverse.eps = 1e-08,max.iterations = 500){
  if(lambda1[1] < 0 | lambda1[1] > 1)
    return(rep(NA,length(x)))
  if(mixture)
    return((1 - lambda1[1]) * (prob * dgl(x = x,lambda1 = lambda1[-1],param = param,inverse.eps = inverse.eps,
                                      max.iterations = max.iterations) + (1 - prob) *
                                dgl(x = x,lambda1 = lambda2,param = param,inverse.eps = inverse.eps,
                                    max.iterations = max.iterations)))
  else
    return((1 - lambda1[1]) * dgl(x = x,lambda1 = lambda1[-1],param = param,inverse.eps = inverse.eps,
                                 max.iterations = max.iterations))
}
