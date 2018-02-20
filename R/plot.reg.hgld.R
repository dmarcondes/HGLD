#' @import GLDEX
#' @import ggplot2
#' @import grDevices
#' @export
#' @title Predict density plot of a Hurdle Generalized Lambda Regression
#'
#' @description Predict density plot of a Hurdle Generalized Lambda Regression.
#'
#' @details Given an object of class \link[HGLD]{reg.hgld} and new values for the covariates,
#'  returns the density function of the fitted HGLD, given the covariates. The contrast on the
#'  \emph{newvalues} data frame must be same of the \emph{data} used in the \link[HGLD]{reg.hgld} object. All
#'  density curves are plotted in the same plot.
#'
#' @param x An object of class \link[HGLD]{reg.hgld}.
#' @param newvalues A data frame with the new values of the covariates. Column names must match the ones given in formulas \emph{loc.formula} and \emph{zero.formula}.
#' @param name A vector with the names of the new values rows. Default is the row names of the new values data frame.
#' @param color The color of each density curve. Must have one color for each row of new values.
#' @param xlab Label of the x-axis.
#' @param xlim A vector with the limits of the x-axis.
#' @param title Legend title.
#' @param ... Arguments to be passed to methods.
#' @return \item{plot}{\link[ggplot2]{ggplot} density plot for the given new values}
#' @examples
#' set.seed(10)
#' tmp <- na.omit(healthcare)
#' data <- tmp[sample(1:nrow(tmp),100),]
#' formula <- log_expense ~ age + sex + log_previous_expense
#' reg <- suppressWarnings(reg.hgld(data,formula,formula,TRUE,n.simu = 10,param = "rs",plotKS = TRUE))
#' newvalues <- tmp[sample(1:nrow(tmp),5),c(2,3,5,8)]
#' plot(reg,newvalues)
#'
#' @references Marcondes, D.; Peixoto, C.; Maia, A. C.; Fitting a Hurdle Generalized Lambda Distribution to healthcare expenses. (2017) \emph{arxiv1712.02183}

plot.reg.hgld <- function(x,newvalues,name = row.names(newvalues),color = NULL,xlab = "Data",
                           xlim = c(min(x$NZdata[[all.vars(x$loc.formula)[1]]]),max(x$NZdata[[all.vars(x$loc.formula)[1]]])),
                           title = "",...){
  # Graphical parameters
  titles <- theme(strip.text = element_text(size = 12), axis.text = element_text(size = 10,
                  color = "black"), axis.title = element_text(size = 12), legend.text = element_text(size = 10),
                  legend.title = element_text(size = 12), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), panel.border = element_blank(),
                  panel.background = element_rect(fill="white",size=0.5, linetype="solid",color = "black"),
                  #legend.position = c(0.85, 0.9),
                  legend.background = element_rect(fill="white",size=0.5, linetype="solid",color = "black"))
  themes <- theme_linedraw()

  # Colors
  if(is.null(color))
    color <- colorRampPalette(c("blue", "red"))(nrow(newvalues))

  # Estimated coefficients
  reg <- x
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

  # Predictors
  newvalues[[all.vars(reg$loc.formula)[1]]] <- rep(1,nrow(newvalues))
  W <- model.matrix(lm(reg$loc.formula,data = newvalues))
  Z <- model.matrix(lm(reg$zero.formula,data = newvalues))
  local <- W %*% beta
  odds <- exp(Z %*% gama)
  lambda0 <- odds/(1+odds)

  # Plotting densitities
  p <- {ggplot(data.frame(x = reg$NZdata[[all.vars(reg$loc.formula)[1]]]),aes(x = x)) + themes + titles + ylab("Density") + xlab(xlab) +
    xlim(xlim)}
  for(i in 1:nrow(newvalues)){
    loop_input = paste("stat_function(fun = function(x) dhgld(x = x - ",local[i],",lambda1 = c(",lambda0[i],",lambda),param = reg$param),aes(colour = '",name[i],"'))", sep="")
    p <- p + eval(parse(text=loop_input))
  }
  names(color) <- name
  p <- p + scale_colour_manual(title,values = color,breaks = rownames(newvalues))
  return(list("plot" = p))
}

