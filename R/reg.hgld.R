#' @import GLDreg
#' @import gamlss
#' @import gamlss.dist
#' @import ggplot2
#' @import gridExtra
#' @import grid
#' @export
#' @title Fit a Hurdle Generalized Lambda Distribution Regression model
#'
#' @description Fit a Hurdle Generalized Lambda Distribution Regression model to a dataset.
#'
#' @details Given a dataset, estimate by the Numerical Maximum Likelihood method the regression coefficients of the model and
#' the five parameters of the error GLD. The regression coefficients that model the location of the distribution are estimated by
#' the functions \link[GLDreg]{GLD.lm} and \link[GLDreg]{GLD.lm.full} of package \link[GLDreg]{GLDreg}. The regression coefficients
#' that model the hurdle parameter of the distribution are estimated by function \link[gamlss]{gamlss}.
#'
#' @param data A dataset containing the variables of the model.
#' @param zero.formula A symbolic expression of the model to be fitted to the hurdle parameter of the distribution.
#' @param loc.formula A symbolic expression of the model to be fitted to the location of the distribution.
#' @param full Whether a simulation method must be applied to derive a confidence interval for the location regression coefficients.
#' @param param "fmkl" or "rs".
#' @param maxit Maximum number of iterations for numerical optimization.
#' @param init Choose a different set of initial values to start the optimization process. This can either be full set of parameters including
#' GLD parameter estimates, or it can just be the coefficient estimates of the regression model.
#' @param alpha Significance level of the Confidence Interval for the GLD regression.
#' @param n.simu Number of times to repeat the simulation runs, defaults to 1000.
#' @param plotKS Whether to plot the KS resample test result within each plot.
#' @param h.bins Number of bins for the GLD Regression normalized quantiles residuals histogram.
#' @return \item{coefficients}{The estimated coefficients of the HGLD regression.}
#' @return \item{Zplot}{A function that generates the four diagnostic plots of the logistic regression.}
#' @return \item{Zres.sumarry}{Summary of the normalised quantile residuals of the logistic regression.}
#' @return \item{Zfit}{Normalised quantile residuals versus fitted values for the logistic regression.}
#' @return \item{Zindex}{Normalised quantile residuals versus index for the logistic regression.}
#' @return \item{Zdensity}{Normalised quantile residuals density for the logistic regression.}
#' @return \item{Zqq}{QQ-norm plot of the normalised quantile residuals for the logistic regression.}
#' @return \item{NZqq}{QQ-plot for the GLD regression residuals.}
#' @return \item{NZquant}{Quantile plot for the GLD regression residuals.}
#' @return \item{NZhistogram}{The histogram of the GLD regression residuals.}
#' @return \item{NZplot}{The normalised quantile residuals plots for the GLD regression.}
#' @return \item{NZres.sumarry}{Summary of the normalised quantile residuals of the GLD regression.}
#' @return \item{NZfit}{Normalised quantile residuals versus fitted values for the GLD regression.}
#' @return \item{NZindex}{Normalised quantile residuals versus index for the GLD regression.}
#' @return \item{NZdensity}{normalised quantile residuals density for the GLD regression.}
#' @return \item{NZqqQuant}{QQ-norm plot of the normalised quantile residuals for the GLD regression.}
#' @return \item{KS}{KS test p-value for the GLD regression.}
#' @return \item{gamlss}{The \link[gamlss]{gamlss} object of the fitted logistic regression.}
#' @return \item{GLDreg}{The \link[GLDreg]{GLDreg} object of the fitted GLD regression.}
#' @return \item{Zdata}{The data used to fit the logistic regression}
#' @return \item{NZdata}{The data used to fit the GLD regression.}
#' @return \item{zero.formula}{A symbolic expression of the model to be fitted to the hurdle parameter of the distribution.}
#' @return \item{loc.formula}{A symbolic expression of the model to be fitted to the location of the distribution.}
#' @return \item{param}{"fmkl" or "rs".}
#' @return \item{full}{Whether a simulation method must be applied to derive a confidence interval for the location regression coefficients.}
#' @examples
#' set.seed(100)
#' tmp <- na.omit(healthcare)
#' data <- tmp[sample(1:nrow(tmp),50),]
#' formula <- log_expense ~ age + sex + log_previous_expense
#' reg <- suppressWarnings(reg.hgld(data = data,zero.formula = formula,loc.formula = formula,
#'                                  full = FALSE,param = "rs"))
#'
#' @references Marcondes, D.; Peixoto, C.; Maia, A. C.; A Survey of a Hurdle Model for Heavy-Tailed Data Based on the Generalized Lambda Distribution. (2017) \emph{arxiv1712.02183}
#' @references Su, S.; Flexible Parametric Quantile Regression Model. (2015), Statistics & Computing May 2015, Volume 25, Issue 3, pp 635-650

reg.hgld <- function(data,zero.formula,loc.formula,full = FALSE,param = "fmkl",maxit = 20000,init = NULL,
                      alpha = 0.05,n.simu = 1000,plotKS = TRUE,h.bins = 50){
  # No visible binding for global variable correction
  x = y = NULL

  # Graphical parameters
  titles <- theme(strip.text = element_text(size = 12), axis.text = element_text(size = 10,
                  color = "black"), axis.title = element_text(size = 12), legend.text = element_text(size = 10),
                  legend.title = element_text(size = 12), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), panel.border = element_blank(),
                  panel.background = element_rect(fill="white",size=0.5, linetype="solid",color = "black"))
  themes <- theme_linedraw()


  # Initializing data
  varsL <- all.vars(loc.formula)
  varsZ <- all.vars(zero.formula)
  Zdata <- data
  Zdata[[varsZ[1]]] <- ifelse(Zdata[[varsZ[1]]] == 0,1,0)
  Ndata <- data[data[[varsZ[1]]] != 0,]

  # Fitting a logistic regression
  Z <- gamlss(formula = zero.formula,family = BI,data = Zdata)
  s <- summary(Z)

  # Fitting a GLD regression
  if(full & param == "rs"){
    NZ <- GLD.lm.full(formula = loc.formula,data = Ndata,param = param,maxit = maxit,fun = fun.RPRS.ml.m,method = "Nelder-Mead",
                      n.simu = n.simu, init = init, summary.plot = FALSE)
  }
  if(full & param == "fmkl"){
    NZ <- GLD.lm.full(formula = loc.formula,data = Ndata,param = param,maxit = maxit,fun = fun.RMFMKL.ml.m,method = "Nelder-Mead",
                      n.simu = n.simu, init = init, summary.plot = FALSE)
  }
  if(!full & param == "rs"){
    NZ <- GLD.lm(formula = loc.formula,data = Ndata,param = param,maxit = maxit,fun = fun.RPRS.ml.m,method = "Nelder-Mead",
                 diagnostics = FALSE,alpha = 0.05, init = init)
  }
  if(!full & param == "fmkl"){
    NZ <- GLD.lm(formula = loc.formula,data = Ndata,param = param,maxit = maxit,fun = fun.RMFMKL.ml.m,method = "Nelder-Mead",
                 diagnostics = FALSE,alpha = 0.05, init = init)
  }

  # Inference about the parameters
  if(full){
    coef <- names(data.frame(NZ[[3]]))
    estimated <- NZ[[1]]$`Estimated parameters`
    k <- data.frame(quantile(NZ[[3]][,1],probs = c(alpha/2),type = 8),quantile(NZ[[3]][,1],probs = 1-alpha/2,type = 8))

    for(i in 2:length(coef))
      k <- rbind.data.frame(k,c(quantile(NZ[[3]][,i],probs = alpha/2,type = 8),quantile(NZ[[3]][,i],probs = 1-alpha/2,type = 8)))
    k <- data.frame(estimated[1:length(coef)],k)
    colnames(k) <- c("Estimate",paste("Lower Bound (",100 * alpha/2,"%)",sep = ""),paste("Upper Bound (",100 *(1 - alpha/2),"%)",sep = ""))
    row.names(k) <- coef
  }
  else{
    estimated <- NZ$`Estimated parameters`
    estimated <- estimated[1:(length(estimated) - 4)]
    coef <- names(estimated)
    k <- data.frame(estimated,NA,NA)
    colnames(k) <- c("Estimate",paste("Lower Bound (",100 * alpha/2,"%)",sep = ""),paste("Upper Bound (",100 *(1 - alpha/2),"%)",sep = ""))
    row.names(k) <- coef
  }
  s <- list("Logistic" = s,"GLD" = k)

  # Diagnostics for the logistic regression
  pfit <- {ggplot(data.frame("y" = Z$residuals,"x" = fitted(Z)),aes(x = x,y = y)) + geom_point() + themes + titles +
      xlab("Fitted Values") + ylab("Quantile Residuals") + ggtitle("Against Fitted values")}
  pindex <- {ggplot(data.frame("y" = Z$residuals,"x" = c(1:length(Z$residuals))),aes(x = x,y = y)) + geom_point() + themes +
      titles + xlab("Index") + ylab("Quantile Residuals") + ggtitle("Against Index")}
  pdensity <- {ggplot(data.frame("x" = Z$residuals),aes(x = x)) + geom_density() +
      geom_point(aes(x = x,y = 0),inherit.aes = FALSE,color = "red") + themes + titles + ylab("Density") +
      xlab("Quantile Residuals") + ggtitle("Density Estimate")}
  pqq <- {ggplot(data.frame("x" = Z$residuals),aes(sample = x)) + stat_qq() + themes + titles + xlab("Theoretical Quantiles") +
      ylab("Sample Quantiles") + geom_abline(intercept = 0, slope = 1, alpha = 0.5) + ggtitle("Normal QQ-Plot")}
  Zplot <- grid.arrange(pfit,pindex,pdensity,pqq)

  qq <- as.data.frame(qqnorm(Z$residuals, plot = FALSE))
  Filliben <- cor(qq$y,qq$x)
  Zres <- data.frame("Value" = c(mean(Z$residuals),var(Z$residuals),skw(Z$residuals),kts(Z$residuals) + 3,
                       Filliben))
  row.names(Zres) <- c("mean","variance","coef. of skewness","coef. of kurtosis","Filliben correlation coefficient")

  # Diagnostic for the GLD
  Ndiag <- diag.reghgld(reg = NZ,param = param,plotKS = plotKS,full = full,h.bins = h.bins)

  # Return
  r <- list("coefficients" = s,"Zplot" = function(){grid.draw(Zplot)},"Zres.summary" = Zres,"Zfit" = pfit,"Zindex" = pindex,
            "Zdensity" = pdensity,"Zqq" = pqq,"NZqq" = Ndiag$qq,"NZquant" = Ndiag$quant,"NZhistogram" = Ndiag$histogram,
            "NZplot" = function(){grid.draw(Ndiag$NZplot)},"NZres.summary" = Ndiag$NZres,"NZfit" = Ndiag$NZfit,"NZindex" = Ndiag$NZindex,
            "NZdensity" = Ndiag$NZdensity,"NZqqQuant" = Ndiag$NZqq,
            "KS" = Ndiag$KS,"gamlss" = Z,"GLDreg" = NZ,"Zdata" = Zdata,"NZdata" = Ndata,"loc.formula" = loc.formula,
            "zero.formula" = zero.formula,"param" = param,"full" = full)
  class(r) <- "reg.hgld"
  return(r)

}
