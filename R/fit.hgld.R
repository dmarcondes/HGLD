#' @import GLDEX
#' @import cluster
#' @export
#' @title Fit the Hurdle Generalized Lambda Distribution
#'
#' @description Fit the Hurdle Generalized Lambda Distribution to a dataset by the Numerical Maximum Likelihood Method.
#'
#' @details Given a dataset, estimate by the Numerical Maximum Likelihood Method the five parameters of the HGLD.
#' Fit both the RS and the fmkl parametrizations. Also fit a mixture of HGLDs.
#'
#' @param data A vector of data.
#' @param mixture Whether a mixture of HGLD must be fitted.
#' @param clustering.m Clustering method used in classifying the dataset into two parts when fitting a mixture of HGLD.
#' Valid arguments include clara, fanny and pam from the cluster library, or threshold, if the data must be divided
#' by a threshold. Default is clara. Or a logical vector specifying how data should be split. See
#' \link[GLDEX]{fun.auto.bimodal.pml} for more details.
#' @param threshold The threshold to divide the data if \emph{clustering.m = "threshold"}.
#' @param leap1 Scrambling (0,1,2,3) for the sobol sequence for the first distribution fit when fitting a mixture of
#' HGLD. See scrambling/leap argument for \link[GLDEX]{runif.sobol}, \link[GLDEX]{runif.halton} or
#' \link[GLDEX]{QUnif} of the \link[GLDEX]{GLDEX} package.
#' @param leap2 Scrambling (0,1,2,3) for the sobol sequence for the second distribution fit when fitting a mixture of
#' HGLD. See scrambling/leap argument for \link[GLDEX]{runif.sobol}, \link[GLDEX]{runif.halton} or
#' \link[GLDEX]{QUnif} of the \link[GLDEX]{GLDEX} package.
#' @param fun1 A character string of either "runif.sobol" (default), "runif.halton" or "QUnif" for the first
#' distribution fit when fitting a mixture of HGLD. See \link[GLDEX]{fun.auto.bimodal.pml} for more details.
#' @param fun2 A character string of either "runif.sobol" (default), "runif.halton" or "QUnif" for the second
#' distribution fit when fitting a mixture of HGLD. See \link[GLDEX]{fun.auto.bimodal.pml} for more details.
#' @param rs.leap Scrambling (0,1,2,3) for the sobol sequence for the RS generalized lambda distribution fit. See scrambling/leap
#' argument for \link[GLDEX]{runif.sobol}, \link[GLDEX]{runif.halton} or \link[GLDEX]{QUnif} of the
#' \link[GLDEX]{GLDEX} package. See \link[GLDEX]{fun.data.fit.ml} for more details.
#' @param fmkl.leap Scrambling (0,1,2,3) for the sobol sequence for the fmkl generalized lambda distribution fit. See scrambling/leap
#' argument for \link[GLDEX]{runif.sobol}, \link[GLDEX]{runif.halton} or \link[GLDEX]{QUnif} of the
#' \link[GLDEX]{GLDEX} package. See \link[GLDEX]{fun.data.fit.ml} for more details.
#' @param rs.init Initial values (lambda3 and lambda4) for the RS generalized lambda distribution. See
#' \link[GLDEX]{fun.data.fit.ml} for more details.
#' @param fmkl.init Initial values (lambda3 and lambda4) for the fmkl generalized lambda distribution. See
#' \link[GLDEX]{fun.data.fit.ml} for more details.
#' @param FUN A character string of either "runif.sobol" (default), "runif.halton" or "QUnif". See
#' \link[GLDEX]{fun.data.fit.ml} for more details.
#' @param no Number of initial random values to find the best initial values for optimization. See
#' \link[GLDEX]{fun.data.fit.ml} for more details.
#' @return \item{par}{The estimate of the HGLD five parameters for both the RS and fmkl parametrizations if \emph{mixture
#' = FALSE}. Otherwise present the estimation of ten parameters: the zero probability mass, the four parameters of
#' each GLD and the clustering parameter p.}
#' @return \item{data}{The data used in the fit.}
#' @return \item{mixture}{Whether a mixture of HGLD was fitted.}
#' @examples
#' set.seed(100)
#' data <- healthcare[sample(1:nrow(healthcare),30),]
#' fit <- fit.hgld(data$log_expense)
#'
#' #mixture
#' set.seed(100)
#' data <- c(rcauchy(20,location = 10),rep(0,10),rcauchy(20))
#' fit2 <- fit.hgld(data = data,mixture = TRUE)
#'
#' @references Marcondes, D.; Peixoto, C.; Maia, A. C.; A Survey of a Hurdle Model for Heavy-Tailed Data Based on the Generalized Lambda Distribution. (2017) \emph{arxiv1712.02183}
#' @references Su, S.; Fitting Single and Mixture of Generalized Lambda Distributions to Data via Discretized and Maximum Likelihood Methods: GLDEX in R.(2007), Journal of Statistical Software: *21* 9.

fit.hgld <- function(data,mixture = FALSE,clustering.m = clara,threshold = NULL,leap1 = 3,leap2 = 3,fun1 = "runif.sobol",
                     fun2 = "runif.sobol",rs.leap = 3,fmkl.leap = 3,rs.init = c(-1.5, 1.5),fmkl.init = c(-0.25, 1.5),
                     FUN = "runif.sobol",no = 10000){
  # Data to be used
  lambda0 <- sum(data == 0)/length(data)
  nzero <- data[data != 0]

  # Fitting the GLD
  if(mixture){
    if(!identical(clustering.m,"threshold")){
      fit.RS <- fun.auto.bimodal.pml(data = nzero,init1.sel = "rprs",init2.sel = "rprs",clustering.m = clustering.m,
                                     init1 = rs.init,init2 = rs.init,leap1 = leap1,leap2 = leap2,fun1 = fun1,
                                     fun2 = fun2,no = no)
      fit.fmkl <- fun.auto.bimodal.pml(data = nzero,init1.sel = "rmfmkl",init2.sel = "rmfmkl",clustering.m = clustering.m,
                                       init1 = fmkl.init,init2 = fmkl.init,leap1 = leap1,leap2 = leap2,fun1 = fun1,
                                       fun2 = fun2,no = no)
    }
    else{
      d1 <- nzero[nzero <= threshold]
      d2 <- nzero[nzero >= threshold]
      p <- length(d1)/(length(d1) + length(d2))
      fit.1 <- fun.data.fit.ml(d1,rs.leap = rs.leap, fmkl.leap = fmkl.leap, rs.init = rs.init,
                               fmkl.init = fmkl.init,FUN = FUN,no = no)[,1:2]
      fit.2 <- fun.data.fit.ml(d2,rs.leap = rs.leap, fmkl.leap = fmkl.leap, rs.init = rs.init,
                               fmkl.init = fmkl.init,FUN = FUN,no = no)[,1:2]
    }
  }
  else
  lambda <- fun.data.fit.ml(nzero,rs.leap = rs.leap, fmkl.leap = fmkl.leap, rs.init = rs.init,
                            fmkl.init = fmkl.init,FUN = FUN,no = no)[,1:2]

  # Returning
  if(mixture){
    if(!identical(clustering.m,"threshold"))
      par <- data.frame("par" = c("lambda0","f1.lambda1","f1.lambda2","f1.lambda3","f1.lambda4","f2.lambda1","f2.lambda2",
                       "f2.lambda3","f2.lambda4","p"),"RS" = c(lambda0,fit.RS$par),
                       "fmkl" = c(lambda0,fit.fmkl$par))
    else
      par <- data.frame("par" = c("lambda0","f1.lambda1","f1.lambda2","f1.lambda3","f1.lambda4","f2.lambda1","f2.lambda2",
                                  "f2.lambda3","f2.lambda4","p"),"RS" = c(lambda0,fit.1[,1],fit.2[,1],p),
                        "fmkl" = c(lambda0,fit.1[,2],fit.2[,2],p))
  }

  else
    par <- data.frame("par" = c("lambda0","lambda1","lambda2","lambda3","lambda4"),"RS" = c(lambda0,lambda[,1]),
                      "fmkl" = c(lambda0,lambda[,2]))
  print(par)
  fit <- list("par" = par,"data" = data,"mixture" = mixture)
  class(fit) <- "fit.hgld"
  return(fit)
}

