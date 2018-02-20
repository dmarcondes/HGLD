#' @import GLDEX
#' @import ggplot2
#' @export
#' @title Plot the Hurdle Generalized Lambda Distribution fit
#'
#' @description Density plot of the fit given by an object of class \link[HGLD]{fit.hgld}.
#'
#' @details The density may be plotted by itself or superimposed by the data histogram. Plot only the nonzero data values.
#'
#' @param x A \link[HGLD]{fit.hgld} object.
#' @param histogram Logical. Whether the estimated density must be superimposed by the data histogram.
#' @param hcolor Color of the histogram.
#' @param hfill Color to fill the histogram.
#' @param bins Number of histogram bins.
#' @param dcolor Color of the RS and FMKL density plot.
#' @param dtype Type of the RS and FMKL density line.
#' @param xlab Label of the x-axis.
#' @param ylab Label of the y-axis.
#' @param ... Arguments to be passed to methods.
#' @return \item{RS}{\link[ggplot2]{ggplot} plot of the fitted RS GLD.}
#' @return \item{FMKL}{\link[ggplot2]{ggplot} plot of the fitted FMKL GLD.}
#' @return \item{RF}{\link[ggplot2]{ggplot} plot of the fitted RS and FMKL GLD on the same plot.}
#' @examples
#' set.seed(100)
#' data <- healthcare[sample(1:nrow(healthcare),100),]
#' fit <- fit.hgld(data$log_expense)
#' plot(fit)
#'
#' #mixture
#' set.seed(100)
#' data <- c(rcauchy(50,location = 10),rep(0,30),rcauchy(50))
#' fit <- fit.hgld(data = data,mixture = TRUE)
#' plot(fit)
#'
#' @references Marcondes, D.; Peixoto, C.; Maia, A. C.; Fitting a Hurdle Generalized Lambda Distribution to healthcare expenses. (2017) \emph{arxiv1712.02183}

plot.fit.hgld <- function(x,histogram = TRUE,hcolor = "white",hfill = "black",bins = 50,dcolor = c("black","red"),
                          dtype = c("solid","solid"),xlab = "Data",ylab = "Density",...){
  # No visible biding for variable correction
  ..density.. <- NULL

  # Graphical parameters
  titles <- theme(strip.text = element_text(size = 12), axis.text = element_text(size = 10,
                  color = "black"), axis.title = element_text(size = 12), legend.text = element_text(size = 10),
                  legend.title = element_text(size = 12), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), panel.border = element_blank(),
                  panel.background = element_rect(fill="white",size=0.5, linetype="solid",color = "black"),
                  legend.position = c(1,1),legend.justification = c(1,1),
                  legend.background = element_rect(fill="white",size=0.5, linetype="solid",color = "black"))
  themes <- theme_linedraw()

  # Data to be used and raw plot
  data <- data.frame(x = x$data[x$data != 0])
  plot <- ggplot(data,aes(x = x))
  mixture <- x$mixture
  lambdaRS1 <- x$par[1:5,2]
  lambdaFMKL1 <- x$par[1:5,3]

  if(mixture){
    pRS <- x$par[10,2]
    lambdaFMKL2 <- x$par[6:9,3]
    lambdaRS2 <- x$par[6:9,2]
    pFMKL <- x$par[10,3]
  }
  else{
    pRS <- NULL
    lambdaFMKL2 <- NULL
    lambdaRS2 <- NULL
    pFMKL <- NULL
  }

  # Densities to plot
  dRS <- function(x){
    dhgld(x = x,mixture = mixture,lambda1 = lambdaRS1,lambda2 = lambdaRS2,prob = pRS,param = "rs")
  }
  dFMKL <- function(x){
    dhgld(x = x,mixture = mixture,lambda1 = lambdaFMKL1,lambda2 = lambdaFMKL2,prob = pFMKL,param = "fmkl")
  }

  # Legend items
  if(mixture){
    nFMKL <- paste("FMKL(",round(lambdaFMKL1[1],2),",",round(lambdaFMKL1[2],2),",",round(lambdaFMKL1[3],2),",",
                   round(lambdaFMKL1[4],2),",",round(lambdaFMKL1[5],2),") \n","FMKL(",round(lambdaFMKL1[1],2),",",round(lambdaFMKL2[1],2),",",
                   round(lambdaFMKL2[2],2),",",round(lambdaFMKL2[3],2),",",round(lambdaFMKL2[4],2),") \np = ",round(pFMKL,2),"\n",sep = "")
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

  # Plots
  if(histogram){
    p1 <- (1 - lambdaFMKL1[1])
    plotRS <- {plot + geom_histogram(aes_string(x = "x",y = paste("..density.. *",p1)),color = hcolor,fill = hfill,bins = bins) +
        stat_function(fun = dRS,linetype = dtype[1],aes(colour = nRS)) + xlab(xlab) + ylab("Density") + themes + titles +
        scale_colour_manual("",values = dcolor[1]) + theme(legend.title = element_blank())}

    plotFMKL <- {plot + geom_histogram(aes_string(x = "x",y = paste("..density.. *",p1)),color = hcolor,fill = hfill,bins = bins) +
        stat_function(fun = dFMKL,linetype = dtype[2],aes(colour = nFMKL)) + xlab(xlab) + ylab("Density") + themes + titles +
        scale_colour_manual("",values = dcolor[2]) + theme(legend.title = element_blank())}

    plotBOTH <- {plot + geom_histogram(aes_string(x = "x",y = paste("..density.. *",p1)),color = hcolor,fill = hfill,bins = bins) +
      stat_function(fun = dRS,linetype = dtype[1],aes(colour = nRS)) +
      stat_function(fun = dFMKL,linetype = dtype[2],aes(colour = nFMKL)) +
      xlab(xlab) + ylab("Density") + themes + titles +
      scale_colour_manual("",values = c(dcolor[2],dcolor[1])) + theme(legend.title = element_blank())}
  }

  else{
    plotRS <- {plot +
        stat_function(fun = dRS,linetype = dtype[1],aes(colour = nRS)) + xlab(xlab) + ylab("Density") + themes + titles +
        scale_colour_manual("",values = dcolor[1]) + theme(legend.title = element_blank())}

    plotFMKL <- {plot +
        stat_function(fun = dFMKL,linetype = dtype[2],aes(colour = nFMKL)) + xlab(xlab) + ylab("Density") + themes + titles +
        scale_colour_manual("",values = dcolor[2]) + theme(legend.title = element_blank())}

    plotBOTH <- {plot +
        stat_function(fun = dRS,linetype = dtype[1],aes(colour = nRS)) +
        stat_function(fun = dFMKL,linetype = dtype[2],aes(colour = nFMKL)) +
        xlab(xlab) + ylab("Density") + themes + titles +
        scale_colour_manual("",values = c(dcolor[2],dcolor[1])) + theme(legend.title = element_blank())}
  }
  return(list("RS" = plotRS,"FMKL" = plotFMKL,"RF" = plotBOTH))
}
