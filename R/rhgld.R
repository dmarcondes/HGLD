#' @import GLDEX
#' @import gamlss
#' @import stats
#' @export
#' @title The Hurdle Generalized Lambda Distribution Family
#'
#' @description Sample from a Hurdle Generalized Lambda Distribution.
#'
#' @details If the parametrization of the RS or FMKL GLD is not acceptable, the function returns NA. This function is based on the \link[GLDEX]{GLDEX} package.
#'
#' @param n	Number of observations to be generated.
#' @param mixture Whether to give the density of a mixture HGLD.
#' @param lambda1 A vector of length 5 with the five parameters of the HGLD or the first HGLD if mixture = TRUE.
#' @param lambda2 A vector of length 4 with the four parameters of the second HGLD.
#' @param prob The cluster parameter for the mixture HGLD.
#' @param param "fmkl" or "rs".
#' @return A sample of a HGLD.
#' @examples
#' library(ggplot2)
#' set.seed(100)
#' qplot(rhgld(n = 1000,lambda1 = c(0.05,0,1,3,6),param = "fmkl"),..density..,geom = "histogram",
#'       bins = 30)
#' qplot(rhgld(n = 1000,mixture = TRUE,lambda1 = c(0.05,0,1,3,6),lambda2 = c(0.3,2,3,6),prob = 0.5,
#'             param = "fmkl"),geom = "histogram",bins = 50)
#'
#' @references Marcondes, D.; Peixoto, C.; Maia, A. C.; A Survey of a Hurdle Model for Heavy-Tailed Data Based on the Generalized Lambda Distribution. (2017) \emph{arxiv1712.02183}
#' @references Su, S.; Fitting Single and Mixture of Generalized Lambda Distributions to Data via Discretized and Maximum Likelihood Methods: GLDEX in R. (2007), Journal of Statistical Software: *21* 9.

rhgld <- function(n,mixture = FALSE,lambda1,lambda2 = NULL,prob = NULL,param = "fmkl"){
  if(mixture){
    z <- 1 - rbinom(n = n,size = 1,prob = lambda1[1])
    zp <- rbinom(n = n,size = 1,prob = prob)
    s1 <- rgl(n = n,lambda1 = lambda1[-1],param = param)
    s2 <- rgl(n = n,lambda1 = lambda2,param = param)
    s <- z * (zp * s1 + (1-zp) * s2)
    return(s)
  }
  else{
    z <- 1 - rbinom(n = n,size = 1,prob = lambda1[1])
    s <- rgl(n = n,lambda1 = lambda1[-1],param = param)
    s <- s*z
  return(s)
  }
}
