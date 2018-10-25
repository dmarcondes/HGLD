#' @import GLDEX
#' @import stats
#' @import utils
#' @export
#' @title The Hurdle Generalized Lambda Distribution Family
#'
#' @description Quantile function of a Hurdle Generalized Lambda Distribution.
#'
#' @details The HGLD must be non-negative, otherwise the quantile is only approximate. If the parametrization of
#' the RS or FMKL GLD is not acceptable, the function returns NA. This function is based on the \link[GLDEX]{GLDEX} package.
#'
#' @param p Vector of probabilities.
#' @param mixture Whether to give the quantile function of a mixture of HGLDs.
#' @param lambda1 A vector of length 5 with the five parameters of the HGLD, or of the first HGLD if mixture = TRUE.
#' @param lambda2 A vector of length 4 with the four parameters of the second HGLD if mixture = TRUE.
#' @param prob The cluster parameter for the mixture HGLD.
#' @param param "fmkl" or "rs".
#' @param trace Whether a progress bar must be printed in order to trace the algorithm.
#' @param inverse.eps Accuracy of calculation for the numerical determination of F(x), defaults to 1e-8.
#' @param max.iterations Maximum number of iterations in the numerical determination of F(x), defaults to 500.
#' @return The quantile function of a HGLD.
#' @examples
#' qhgld(p = seq(0.05,1,0.05),lambda1 = c(0.540,3.561,0.019,0.009,0.022),param = "rs")
#'
#' #mixture
#' qhgld(p = seq(0.05,1,0.05),mixture = TRUE,lambda1 = c(0.1,8,1,3,6),lambda2 = c(0.3,10,3,6),
#'       prob = 0.5,param = "fmkl")
#'
#' @references Marcondes, D.; Peixoto, C.; Maia, A. C.; A Survey of a Hurdle Model for Heavy-Tailed Data Based on the Generalized Lambda Distribution. (2017) \emph{arxiv1712.02183}
#' @references Su, S.; Fitting Single and Mixture of Generalized Lambda Distributions to Data via Discretized and Maximum Likelihood Methods: GLDEX in R. (2007), Journal of Statistical Software: *21* 9.

qhgld <- function(p,mixture = FALSE,lambda1,lambda2 = NULL,prob = NULL,param = "fmkl",trace = FALSE,
                  inverse.eps = 1e-08,max.iterations = 500){
  CDF <- function(q) phgld(q = q,mixture = mixture,lambda1 = lambda1,lambda2 = lambda2,prob = prob,param = param,
                           inverse.eps = inverse.eps,max.iterations = max.iterations)
  r <- range.hgld(mixture = mixture,lambda1 = lambda1,lambda2 = lambda2,param = param)
  l <- ifelse(is.finite(r[1]),r[1],-1e106)
  u <- ifelse(is.finite(r[2]),r[2],1e106)
  p <- data.frame(p = p,order = 1:length(p))
  p <- p[order(p$p),]

  quant <- function(v,lower,upper){
    uniroot(function(x) (CDF(x) - v), lower = lower, upper = upper)[1]
  }

  if(trace)
    pb <- txtProgressBar(min = 0, max = nrow(p), style = 3)
  quantile <- vector()
  for(i in 1:nrow(p)){
    v <- p$p[i]
    if(v == 0)
      quantile[i] <- r[1]
    else if(v == 1)
      quantile[i] <- r[2]
    else if(i == 1 & v != 0 & v != 1)
      quantile[i] <- quant(v,lower = l,upper = u)
    else if(!is.finite(unlist(quantile[i-1])))
      quantile[i] <- quant(v,lower = l,upper = u)
    else
      quantile[i] <- quant(v,lower = unlist(quantile[i-1])-0.1,upper = u)

    if(trace)
      setTxtProgressBar(pb,i)
  }
  if(trace)
    cat("\n")
  quantile <- round(unlist(quantile),3)
  names(quantile) <- as.character(p$p)
  quantile <- quantile[order(p$order)]
  return(quantile)
}
