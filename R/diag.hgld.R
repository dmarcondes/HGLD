#' @import GLDEX
#' @import ggplot2
#' @export
#' @title Diagnostic for the Hurdle Generalized Lambda Distribution fit
#'
#' @description Diagnostic plots and measures of goodness-of-fit for a Hurdle Generalized Lambda Distribution fit.
#'
#' @details The diagnostics techniques are applied to the nonzero data values. Returns the qq-plot and the quantile plot between the
#' data and the theoretical fitted HGLD. Also returns a table comparing sample moments with theoretical moments. A Kolmogorov-Simornov resample
#' test is performed and the percentage of the times that the null hypotheses, i.e., goodness-of-fit, is not rejected is displayed.
#' All diagnostics are performed for both the RS and FMKL GLD.
#'
#' @param fit An object of class \link[HGLD]{fit.hgld}.
#' @param facet Whether the plots must be faceted for better visualization.
#' @param facet.breaks The breaks in which to facet the data. Must be the endpoints of the intervals.
#' @param facet.labels The labels of the categories given by the facet breaks.
#' @param facet.ncol The number of columns for the facet plot.
#' @param trace Whether progress bar must be printed in order to trace the algorithm.
#' @param KS Whether the resample KS test must be performed to the nonzero values.
#' @param no.test Total number of KS tests required.
#' @param len Number of data to sample at each KS test.
#' @param alpha Significance level of KS test.
#' @param plotKS Whether to plot the KS resample test result within each plot.
#' @return \item{qqRS}{\link[ggplot2]{ggplot} qq-plot for the RS HGLD.}
#' @return \item{qqFMKL}{\link[ggplot2]{ggplot} qq-plot for the FMKL HGLD.}
#' @return \item{quantRS}{\link[ggplot2]{ggplot} quantile plot for the RS HGLD.}
#' @return \item{quantFMKL}{\link[ggplot2]{ggplot} quantile plot for the FMKL HGLD.}
#' @return \item{moments}{Moments comparison for the GLD fitted to the nonzero data for both parametrizations. Those are not the
#' moments of the HGLD.}
#' @return \item{KS}{Percentage of no rejection for the KS resample test.}
#' @examples
#' set.seed(100)
#' data <- healthcare[sample(1:nrow(healthcare),100),]
#' fit <- fit.hgld(data$log_expense)
#' d <- suppressWarnings(diag.hgld(fit,facet = FALSE,plotKS = FALSE))
#'
#' #mixture
#' set.seed(100)
#' data <- c(rcauchy(50,location = 10),rep(0,30),rcauchy(50))
#' fit <- fit.hgld(data = data,mixture = TRUE)
#' d <- suppressWarnings(diag.hgld(fit))
#'
#' @references Marcondes, D.; Peixoto, C.; Maia, A. C.; A Survey of a Hurdle Model for Heavy-Tailed Data Based on the Generalized Lambda Distribution. (2017) \emph{arxiv1712.02183}
#' @references Su, S. Fitting Single and Mixture of Generalized Lambda Distributions to Data via Discretized and Maximum Likelihood Methods: GLDEX in R. (2007), Journal of Statistical Software: *21* 9.

diag.hgld <- function(fit,facet = FALSE,facet.breaks,facet.labels,facet.ncol,trace = TRUE,no.test = 1000,
                      len = floor(0.9 * length(fit$data[fit$data != 0])),alpha = 0.05,plotKS = TRUE,KS = TRUE){

  # No visible binding for global variable correction
  x = y = NULL

  if(trace)
    cat("Diagnostic for the Hurdle Generalized Lambda Distribution fit \n")

  # Data to be used
  if(KS == FALSE)
    plotKS <- FALSE
  mixture <- fit$mixture
  nonzero <- fit$data[fit$data != 0]
  lambdaRS1 <- fit$par[1:5,2]
  lambdaFMKL1 <- fit$par[1:5,3]

  if(mixture){
    pRS <- fit$par[10,2]
    lambdaFMKL2 <- fit$par[6:9,3]
    lambdaRS2 <- fit$par[6:9,2]
    pFMKL <- fit$par[10,3]
  }
  else{
    pRS <- NULL
    lambdaFMKL2 <- NULL
    lambdaRS2 <- NULL
    pFMKL <- NULL
  }

  if(facet)
    facet.values <- factor(cut(x = fit$data,breaks = facet.breaks,labels = facet.labels,include.lowest = TRUE))
  else
    facet.values <- 1

  # Legend items
  if(mixture){
    nFMKL <- paste("FMKL(",round(lambdaFMKL1[1],2),",",round(lambdaFMKL1[2],2),",",round(lambdaFMKL1[3],2),",",
                   round(lambdaFMKL1[4],2),",",round(lambdaFMKL1[5],2),") \n","FMKL(",round(lambdaFMKL1[1],2),",",round(lambdaFMKL2[1],2),",",
                   round(lambdaFMKL2[2],2),",",round(lambdaFMKL2[3],2),",",round(lambdaFMKL2[4],2),") \np = ",round(pFMKL,2),sep = "")
    nRS <- paste("RS(",round(lambdaRS1[1],2),",",round(lambdaRS1[2],2),",",round(lambdaRS1[3],2),",",
                   round(lambdaRS1[4],2),",",round(lambdaRS1[5],2),") \n","RS(",round(lambdaRS1[1],2),",",round(lambdaRS2[1],2),",",
                   round(lambdaRS2[2],2),",",round(lambdaRS2[3],2),",",round(lambdaRS2[4],2),") \np = ",round(pRS,2),sep = "")
  }
  else{
    nFMKL <- paste("FMKL(",round(lambdaFMKL1[1],2),",",round(lambdaFMKL1[2],2),",",round(lambdaFMKL1[3],2),",",
                   round(lambdaFMKL1[4],2),",",round(lambdaFMKL1[5],2),")",sep = "")
    nRS <- paste("RS(",round(lambdaRS1[1],2),",",round(lambdaRS1[2],2),",",round(lambdaRS1[3],2),",",
                 round(lambdaRS1[4],2),",",round(lambdaRS1[5],2),")",sep = "")
  }

  # Graphical parameters
  titles <- theme(strip.text = element_text(size = 12), axis.text = element_text(size = 10,
                  color = "black"), axis.title = element_text(size = 12), legend.text = element_text(size = 10),
                  legend.title = element_text(size = 12), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), panel.border = element_blank(),
                  panel.background = element_rect(fill="white",size=0.5, linetype="solid",color = "black"),
                  legend.position = c(0,1),legend.justification = c(0, 1),
                  legend.background = element_rect(fill="white",size=0.5, linetype="solid",color = "black"))
  titlesF <- theme(strip.text = element_text(size = 12), axis.text = element_text(size = 10,
                   color = "black"), axis.title = element_text(size = 12), legend.text = element_text(size = 10),
                   legend.title = element_text(size = 12), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), panel.border = element_blank(),
                   panel.background = element_rect(fill="white",size=0.5, linetype="solid",color = "black"),
                   legend.position = "bottom",legend.justification = c(0, 1),
                   legend.background = element_rect(fill="white",size=0.5, linetype="solid",color = "black"))
  themes <- theme_linedraw()

  # QQ-plotRS
  observedRS <- rank(fit$data)/length(fit$data)
  if(trace)
    cat("\n Computing the theoretical quantiles for the RS parametrization \n")
  theoreticalRS <- qhgld(p = observedRS,mixture = mixture,lambda1 = lambdaRS1,lambda2 = lambdaRS2,prob = pRS,
                         param = "rs",trace = trace)
  observedRS <- observedRS[is.finite(theoreticalRS)]
  o <- fit$data[is.finite(theoreticalRS)]
  RSfacet.values <- facet.values[is.finite(theoreticalRS)]
  theoreticalRS <- theoreticalRS[is.finite(theoreticalRS)]
  qqRS <- data.frame(t = theoreticalRS,o = o,facet.values = RSfacet.values)
  if(KS){
    if(trace)
      cat("\n Kolmogorov-Smirnov Tests for the RS parametrization \n")
    ksRS <- ks.test.hgld(data = fit$data,lambda1 = lambdaRS1[-1],lambda2 = lambdaRS2,p = pRS,mixture = mixture,
                          param = "rs",alpha = alpha,len = len,no.test = no.test,trace = trace)
  }

  if(facet)
    plotqqRS <- ggplot(qqRS,aes(x = o,y = t)) + facet_wrap(~ facet.values,ncol = facet.ncol,scales = "free") + themes + titlesF
  else
    plotqqRS <- ggplot(qqRS,aes(x = o,y = t)) + themes + titles

  if(plotKS)
    plotqqRS <- {plotqqRS + geom_point(size = 1) + geom_smooth(color = "black",method = "lm",se = F) +
        xlab("Sample Quantiles") + ylab("Theoretical Quantiles") +
        annotate(geom = "text",x = Inf,y = -Inf,label = paste("K-S resample test = ",ksRS,"%"),hjust = 1, vjust = -1)}
  else
    plotqqRS <- {plotqqRS + geom_point(size = 1) + geom_smooth(color = "black",method = "lm",se = F) +
        xlab("Sample Quantiles") + ylab("Theoretical Quantiles")}

  # Quantile plot RS
  quantRS <- data.frame(x = observedRS,y = o,facet.values = RSfacet.values)

  if(facet)
    plotquantRS <- ggplot(quantRS,aes(x = x,y = y)) + facet_wrap(~ facet.values,ncol = facet.ncol,scales = "free") + themes + titlesF
  else
    plotquantRS <- ggplot(quantRS,aes(x = x,y = y)) + themes + titles

  if(plotKS)
    plotquantRS <- {plotquantRS + geom_line(aes(colour = "Data")) +
        stat_function(fun = function(x) qhgld(p = x,mixture = mixture,lambda1 = lambdaRS1,lambda2 = lambdaRS2,
                                              param = "rs",prob = pRS),aes(colour = nRS)) +
        scale_colour_manual(values = c("black","red")) +
        xlab("Quantile") + ylab("Value") +
        annotate(geom = "text",x = Inf,y = -Inf,label = paste("K-S resample test =",ksRS,"%"),hjust = 1, vjust = -1) +
        theme(legend.title = element_blank())}
  else
    plotquantRS <- {plotquantRS + geom_line(aes(colour = "Data")) +
        stat_function(fun = function(x) qhgld(p = x,mixture = mixture,lambda1 = lambdaRS1,lambda2 = lambdaRS2,
                                              param = "rs",prob = pRS),aes(colour = nRS)) +
        scale_colour_manual(values = c("black","red")) +
        xlab("Quantile") + ylab("Value") + theme(legend.title = element_blank())}

  # QQ-plotFMKL
  observedFMKL <- rank(fit$data)/length(fit$data)
  if(trace)
    cat("\n Computing the theoretical quantiles for the FMKL parametrization \n")
  theoreticalFMKL <- qhgld(p = observedFMKL,mixture = mixture,lambda1 = lambdaFMKL1,lambda2 = lambdaFMKL2,prob = pFMKL,
                           param = "fmkl",trace = trace)
  observedFMKL <- observedFMKL[is.finite(theoreticalFMKL)]
  o <- fit$data[is.finite(theoreticalFMKL)]
  FMKLfacet.values <- facet.values[is.finite(theoreticalFMKL)]
  theoreticalFMKL <- theoreticalFMKL[is.finite(theoreticalFMKL)]
  qqFMKL <- data.frame(t = theoreticalFMKL,o = o,facet.values = FMKLfacet.values)

  if(KS){
    if(trace)
      cat("\n Kolmogorov-Smirnov Tests for the FMKL parametrization \n")
    ksFMKL <- ks.test.hgld(data = fit$data,lambda1 = lambdaFMKL1[-1],lambda2 = lambdaFMKL2,p = pFMKL,mixture = mixture,
                            param = "fmkl",alpha = alpha,len = len,no.test = no.test,trace = trace)
  }


  if(trace)
    cat("\n Just another moment, we are computing the results... \n")

  if(facet)
    plotqqFMKL <- ggplot(qqFMKL,aes(x = o,y = t)) + facet_wrap(~ facet.values,ncol = facet.ncol,scales = "free") + themes + titlesF
  else
    plotqqFMKL <- ggplot(qqFMKL,aes(x = o,y = t)) + themes + titles

  if(plotKS)
    plotqqFMKL <- {plotqqFMKL + geom_point(size = 1) + geom_smooth(color = "black",method = "lm",se = F) +
        xlab("Sample Quantiles") + ylab("Theoretical Quantiles") +
        annotate(geom = "text",x = Inf,y = -Inf,label = paste("K-S resample test =",ksFMKL,"%"),hjust = 1, vjust = -1)}
  else
    plotqqFMKL <- {plotqqFMKL + geom_point(size = 1) + geom_smooth(color = "black",method = "lm",se = F) +
        xlab("Sample Quantiles") + ylab("Theoretical Quantiles")}

  # Quantile plot FMKL
  quantFMKL <- data.frame(x = observedFMKL,y = o,facet.values = FMKLfacet.values)

  if(facet)
    plotquantFMKL <- ggplot(quantFMKL,aes(x = x,y = y)) + facet_wrap(~ facet.values,ncol = facet.ncol,scales = "free") + themes + titlesF
  else
    plotquantFMKL <- ggplot(quantFMKL,aes(x = x,y = y)) + themes + titles

  if(plotKS)
    plotquantFMKL <- {plotquantFMKL + geom_line(aes(colour = "Data")) +
        stat_function(fun = function(x) qhgld(p = x,mixture = mixture,lambda1 = lambdaFMKL1,lambda2 = lambdaFMKL2,
                                              param = "fmkl",prob = pFMKL),aes(colour = nFMKL)) +
        scale_colour_manual(values = c("black","red")) +
        xlab("Quantile") + ylab("Sample") +
        annotate(geom = "text",x = Inf,y = -Inf,label = paste("K-S resample test =",ksFMKL,"%"),hjust = 1, vjust = -1) + theme(legend.title = element_blank())}
  else
    plotquantFMKL <- {plotquantFMKL + geom_line(aes(colour = "Data")) +
        stat_function(fun = function(x) qhgld(p = x,mixture = mixture,lambda1 = lambdaFMKL1,lambda2 = lambdaFMKL2,
                                              param = "fmkl",prob = pFMKL),aes(colour = nFMKL)) +
        scale_colour_manual(values = c("black","red")) +
        xlab("Quantile") + ylab("Sample") + theme(legend.title = element_blank())}

  # Moments comparison
  if(mixture){
    momentsRS <- fun.theo.bi.mv.gld(L1 = lambdaRS1[2],L2 = lambdaRS1[3],L3 = lambdaRS1[4],L4 = lambdaRS1[5],
                                    M1 = lambdaRS2[1],M2 = lambdaRS2[2],M3 = lambdaRS2[3],M4 = lambdaRS2[4],
                                    param1 = "rs",param2 = "rs",p1 = pRS,normalise = "Y")
    momentsFMKL <- fun.theo.bi.mv.gld(L1 = lambdaFMKL1[2],L2 = lambdaFMKL1[3],L3 = lambdaFMKL1[4],L4 = lambdaFMKL1[5],
                                      M1 = lambdaFMKL2[1],M2 = lambdaFMKL2[2],M3 = lambdaFMKL2[3],M4 = lambdaFMKL2[4],
                                      param1 = "fmkl",param2 = "fmkl",p1 = pFMKL,normalise = "Y")
    omomentsRS <- unlist(fun.moments(x = nonzero,normalise = "Y"))
    m <- data.frame("Moments" = c("Mean","Variance","Skewness","Kurtosis"),"Observed" = omomentsRS,"Theoretical RS" = momentsRS,
                    "Theoretical FMKL" = momentsFMKL)
    row.names(m) <- NULL
  }
  else{
    momentsRS <- fun.theo.mv.gld(lambdaRS1[-1],param = "rs",normalise = "Y")
    momentsFMKL <- fun.theo.mv.gld(lambdaFMKL1[-1],param = "fmkl",normalise = "Y")
    omomentsRS <- unlist(fun.moments(x = nonzero,normalise = "Y"))
    m <- data.frame("Moments" = c("Mean","Variance","Skewness","Kurtosis"),"Observed" = omomentsRS,"Theoretical RS" = momentsRS,
                    "Theoretical FMKL" = momentsFMKL)
    row.names(m) <- NULL
  }

  # Return
  if(KS)
    r <- list("qqRS" = plotqqRS,
              "qqFMKL" = plotqqFMKL,
              "quantRS" = plotquantRS,
              "quantFMKL" = plotquantFMKL,
              "moments" = m,
              "KS" = data.frame("RS" = 100*ksRS,"FMKL" = 100*ksFMKL))
  else
    r <- list("qqRS" = plotqqRS,
              "qqFMKL" = plotqqFMKL,
              "quantRS" = plotquantRS,
              "quantFMKL" = plotquantFMKL,
              "moments" = m)
  if(trace)
    cat("\n Finished! \n")
  return(r)
  }
