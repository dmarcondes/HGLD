# Diagnostics for the HGLD regression
diag.reghgld <- function(reg,param,plotKS,full,h.bins){
  x = y = NULL
  ..density.. = NULL
  if(full){
    nonzero <- reg[[1]]$Residual
    lambda <- reg[[1]]$`Estimated parameters`
    lambda <- lambda[(length(lambda)-3):length(lambda)]
    fitted <- reg[[1]]$Fitted
  }
  else{
    nonzero <- reg$Residual
    lambda <- reg$`Estimated parameters`
    lambda <- lambda[(length(lambda)-3):length(lambda)]
    fitted <- reg$Fitted
  }
  nFMKL <- paste("FMKL(",round(lambda[1],2),",",round(lambda[2],2),",",round(lambda[3],2),",",
                 round(lambda[4],2),")",sep = "")
  nRS <- paste("RS(",round(lambda[1],2),",",round(lambda[2],2),",",round(lambda[3],2),",",
               round(lambda[4],2),")",sep = "")

  titles <- theme(strip.text = element_text(size = 12), axis.text = element_text(size = 10,
                  color = "black"), axis.title = element_text(size = 12), legend.text = element_text(size = 10),
                  legend.title = element_text(size = 12), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), panel.border = element_blank(),
                  panel.background = element_rect(fill="white",size=0.5, linetype="solid",color = "black"),
                  legend.position = c(0.15, 0.9),
                  legend.background = element_rect(fill="white",size=0.5, linetype="solid",color = "black"))
  themes <- theme_linedraw()

  if(param == "rs"){
    # qq-plotRS
    observedRS <- rank(nonzero)/length(nonzero)
    theoreticalRS <- qgl(p = observedRS,lambda1 = lambda,param = "rs")
    observedRS <- observedRS[is.finite(theoreticalRS)]
    o <- nonzero[is.finite(theoreticalRS)]
    theoreticalRS <- theoreticalRS[is.finite(theoreticalRS)]
    qqRS <- data.frame("t" = theoreticalRS,"o" = o)
    ks <- round(ks.test(x = nonzero,y = pgl,lambda1 = lambda,param = "rs")$p.value,4)

    if(plotKS)
      plotqq <- {ggplot(qqRS,aes(x = o,y = t)) + geom_point() + geom_smooth(color = "black",method = "lm",se = F) + themes + titles +
          xlab("Residuals Quantiles") + ylab("Theoretical Quantiles") +
          annotate(geom = "text",x = (min(qqRS$o) + 0.1 * abs(max(qqRS$o))),
                   y = max(qqRS$t),label = paste("K-S test p-value = ",ks))}
    else
      plotqq <- {ggplot(qqRS,aes(x = o,y = t)) + geom_point() + geom_smooth(color = "black",method = "lm",se = F) + themes + titles +
          xlab("Residuals Quantiles") + ylab("Theoretical Quantiles")}

    # Quantile plot RS
    quantRS <- data.frame("x" = observedRS,"y" = o)

    if(plotKS)
      plotquant <- {ggplot(quantRS,aes(x = x,y = y)) + geom_line(aes(colour = "Residuals")) +
          stat_function(fun = function(x) qgl(x,lambda1 = lambda,param = "rs"),aes(colour = nRS)) +
          themes + titles + scale_colour_manual(values = c("black","red")) +
          xlab("Quantile") + ylab("Residuals") +
          annotate(geom = "text",x = (min(quantRS$x) + 0.1 * abs(max(quantRS$x))),
                   y = 0.8*max(quantRS$y),label = paste("K-S test p-value =",ks)) + theme(legend.title = element_blank())}
    else
      plotquant <- {ggplot(quantRS,aes(x = x,y = y)) + geom_line(aes(colour = "Residuals")) +
          stat_function(fun = function(x) qgl(x,lambda1 = lambda,param = "rs"),aes(colour = nRS)) +
          themes + titles + scale_colour_manual(values = c("black","red")) +
          xlab("Quantile") + ylab("Residuals") + theme(legend.title = element_blank())}

    # Histogram
    histogram <- ggplot(data.frame("x" = nonzero),aes(x = x)) + themes + titles +
      geom_histogram(color = "black",fill = "white",aes(y = ..density..),bins = h.bins) +
      stat_function(fun = function(x) dgl(x = x,lambda1 = lambda,param = "rs")) + xlab("Residuals") + ylab("Density")

    # Quantile residuals
    quantile_res <- qnorm(p = pgl(q = nonzero,lambda1 = lambda,param = "rs"))
  }

  if(param == "fmkl"){
    # qq-plotFMKL
    observedFMKL <- rank(nonzero)/length(nonzero)
    theoreticalFMKL <- qgl(p = observedFMKL,lambda1 = lambda,param = "fmkl")
    observedFMKL <- observedFMKL[is.finite(theoreticalFMKL)]
    o <- nonzero[is.finite(theoreticalFMKL)]
    theoreticalFMKL <- theoreticalFMKL[is.finite(theoreticalFMKL)]
    qqFMKL <- data.frame("t" = theoreticalFMKL,"o" = o)
    ks <- round(ks.test(x = nonzero,y = pgl,lambda1 = lambda,param = "fmkl")$p.value,4)

    if(plotKS)
      plotqq <- {ggplot(qqFMKL,aes(x = o,y = t)) + geom_point() + geom_smooth(color = "black",method = "lm",se = F) + themes + titles +
          xlab("Residuals Quantiles") + ylab("Theoretical Quantiles") +
          annotate(geom = "text",x = (min(qqFMKL$o) + 0.1 * abs(max(qqFMKL$o))),
                   y = max(qqFMKL$t),label = paste("K-S test p-value =",ks))}
    else
      plotqq <- {ggplot(qqFMKL,aes(x = o,y = t)) + geom_point() + geom_smooth(color = "black",method = "lm",se = F) + themes + titles +
          xlab("Residuals Quantiles") + ylab("Theoretical Quantiles")}

    # Quantile plot FMKL
    quantFMKL <- data.frame("x" = observedFMKL,"y" = o)
    if(plotKS)
      plotquant <- {ggplot(quantFMKL,aes(x = x,y = y)) + geom_line(aes(colour = "Residuals")) +
          stat_function(fun = function(x) qgl(x,lambda1 = lambda,param = "fmkl"),aes(colour = nFMKL)) +
          themes + titles + scale_colour_manual(values = c("red","black")) +
          xlab("Quantile") + ylab("Residuals") +
          annotate(geom = "text",x = (min(quantFMKL$x) + 0.1 * abs(max(quantFMKL$x))),
                   y = 0.8*max(quantFMKL$y),label = paste("K-S test p-value =",ks)) + theme(legend.title = element_blank())}
    else
      plotquant <- plotquant <- {ggplot(quantFMKL,aes(x = x,y = y)) + geom_line(aes(colour = "Residuals")) +
          stat_function(fun = function(x) qgl(x,lambda1 = lambda,param = "fmkl"),aes(colour = nFMKL)) +
          themes + titles + scale_colour_manual(values = c("red","black")) +
          xlab("Quantile") + ylab("Residuals") + theme(legend.title = element_blank())}

    #Histogram
    histogram <- ggplot(data.frame("x" = nonzero),aes(x = x)) + themes + titles +
      geom_histogram(color = "black",fill = "white",aes(y = ..density..),bins = h.bins) +
      stat_function(fun = function(x) dgl(x = x,lambda1 = lambda,param = "fmkl")) + xlab("Residuals") + ylab("Density")

    # Quantile residuals
    quantile_res <- qnorm(p = pgl(q = nonzero,lambda1 = lambda,param = "fmkl"))
  }

  # Quantile residuals plots
  NZfit <- {ggplot(data.frame("y" = quantile_res,"x" = fitted),aes(x = x,y = y)) + geom_point() + themes + titles +
      xlab("Fitted Values") + ylab("Quantile Residuals") + ggtitle("Against Fitted values")}
  NZindex <- {ggplot(data.frame("y" = quantile_res,"x" = c(1:length(quantile_res))),aes(x = x,y = y)) + geom_point() + themes +
      titles + xlab("Index") + ylab("Quantile Residuals") + ggtitle("Against Index")}
  NZdensity <- {ggplot(data.frame("x" = quantile_res),aes(x = x)) + geom_density() +
      geom_point(aes(x = x,y = 0),inherit.aes = FALSE,color = "red") + themes + titles + ylab("Density") +
      xlab("Quantile Residuals") + ggtitle("Density Estimate")}
  NZqq <- {ggplot(data.frame("x" = quantile_res),aes(sample = x)) + stat_qq() + themes + titles + xlab("Theoretical Quantiles") +
      ylab("Sample Quantiles") + geom_abline(intercept = 0, slope = 1, alpha = 0.5) + ggtitle("Normal QQ-Plot")}
  NZplot <- grid.arrange(NZfit,NZindex,NZdensity,NZqq)

  # Quantile residuals statistics
  qq <- as.data.frame(qqnorm(quantile_res, plot = FALSE))
  Filliben <- cor(qq$y,qq$x)
  NZres <- data.frame("Value" = c(mean(quantile_res),var(quantile_res),skw(quantile_res),kts(quantile_res) + 3,
                                 Filliben))
  row.names(NZres) <- c("mean","variance","coef. of skewness","coef. of kurtosis","Filliben correlation coefficient")

  # Return
  r <- list("qq" = plotqq,"NZplot" = NZplot,"NZfit" = NZfit,"NZindex" = NZindex,"NZdensity" = NZdensity,"NZqq" = NZqq,
            "quant" = plotquant,"histogram" = histogram,"NZres" = NZres,
            "KS" = ks)
  return(r)
}

# Skewness function from e1071 package
skw <- function (x, na.rm = FALSE, type = 3){
  if (any(ina <- is.na(x))) {
    if (na.rm)
      x <- x[!ina]
    else return(NA)
  }
  if (!(type %in% (1:3)))
    stop("Invalid 'type' argument.")
  n <- length(x)
  x <- x - mean(x)
  y <- sqrt(n) * sum(x^3)/(sum(x^2)^(3/2))
  if (type == 2) {
    if (n < 3)
      stop("Need at least 3 complete observations.")
    y <- y * sqrt(n * (n - 1))/(n - 2)
  }
  else if (type == 3)
    y <- y * ((1 - 1/n))^(3/2)
  y
}

# Kurtosis function from e1071 package
kts <- function (x, na.rm = FALSE, type = 3){
  if (any(ina <- is.na(x))) {
    if (na.rm)
      x <- x[!ina]
    else return(NA)
  }
  if (!(type %in% (1:3)))
    stop("Invalid 'type' argument.")
  n <- length(x)
  x <- x - mean(x)
  r <- n * sum(x^4)/(sum(x^2)^2)
  y <- if (type == 1)
    r - 3
  else if (type == 2) {
    if (n < 4)
      stop("Need at least 4 complete observations.")
    ((n + 1) * (r - 3) + 6) * (n - 1)/((n - 2) * (n - 3))
  }
  else r * (1 - 1/n)^2 - 3
  y
}

# Modified qhgld
mqhgld <- function(p,lambda,param,local){
  quant <- data.frame()
  for(i in 1:nrow(lambda)){
    l <- as.numeric(as.vector(lambda[i,]))
    q <- vector()
    q[p <= l[1]] <- 0
    if(length(q[p > l[1]]) > 0){
      ptemp <- (p[p > l[1]] - l[1])/(1 - l[1])
      q[p > l[1]] <- qgl(p = ptemp,lambda1 = l[2],lambda2 = l[3],lambda3 = l[4],lambda4 = l[5],param = param)
      q <- unlist(q)
      q[p > l[1]] <- q[p > l[1]] + local[i]
    }
    quant <- rbind.data.frame(quant,q)
  }
  colnames(quant) <- p
  return(quant)
}

# Range of the HGLD
range.hgld <- function(mixture = FALSE,lambda1,lambda2 = NULL,param = "fmkl"){
  lambda0 <- lambda1[1]
  lambda1 <- lambda1[-1]

  if(mixture){
    if(param == "fmkl"){
      if(lambda1[3] > 0)
        infimo1 <- lambda1[1] - (1/(lambda1[2] * lambda1[3]))
      if(lambda1[3] <= 0)
        infimo1 <- -Inf
      if(lambda1[4] > 0)
        sup1 <- lambda1[1] + (1/(lambda1[2] * lambda1[4]))
      if(lambda1[4] <= 0)
        sup1 <- Inf
      if(lambda2[3] > 0)
        infimo2 <- lambda2[1] - (1/(lambda2[2] * lambda2[3]))
      if(lambda2[3] <= 0)
        infimo2 <- -Inf
      if(lambda2[4] > 0)
        sup2 <- lambda2[1] + (1/(lambda2[2] * lambda2[4]))
      if(lambda2[4] <= 0)
        sup2 <- Inf
      infimo <- min(infimo1,infimo2)
      sup <- max(sup1,sup2)
    }

    if(param == "rs"){
      if(min(lambda1[3],lambda1[4]) > 0){
        infimo1 <- lambda1[1] - 1/lambda1[2]
        sup1 <- lambda1[1] + 1/lambda1[2]
      }
      if(lambda1[3] > 0 & lambda1[4] == 0){
        infimo1 <- lambda1[1]
        sup1 <- lambda1[1] + 1/lambda1[2]
      }
      if(lambda1[3] == 0 & lambda1[4] > 0){
        infimo1 <- lambda1[1] - 1/lambda1[2]
        sup1 <- lambda1[1]
      }
      if(max(lambda1[3],lambda1[4]) < 0){
        infimo1 <- -Inf
        sup1 <- Inf
      }
      if(lambda1[3] < 0 & lambda1[4] == 0){
        infimo1 <- -Inf
        sup1 <- lambda1[1] + 1/lambda1[2]
      }
      if(lambda1[3] == 0 & lambda1[4] < 0){
        infimo1 <- lambda1[1] - 1/lambda1[2]
        sup1 <- Inf
      }

      if(min(lambda2[3],lambda2[4]) > 0){
        infimo2 <- lambda2[1] - 1/lambda2[2]
        sup2 <- lambda2[1] + 1/lambda2[2]
      }
      if(lambda2[3] > 0 & lambda2[4] == 0){
        infimo2 <- lambda2[1]
        sup2 <- lambda2[1] + 1/lambda2[2]
      }
      if(lambda2[3] == 0 & lambda2[4] > 0){
        infimo2 <- lambda2[1] - 1/lambda2[2]
        sup2 <- lambda2[1]
      }
      if(max(lambda2[3],lambda2[4]) < 0){
        infimo2 <- -Inf
        sup2 <- Inf
      }
      if(lambda2[3] < 0 & lambda2[4] == 0){
        infimo2 <- -Inf
        sup2 <- lambda2[1] + 1/lambda2[2]
      }
      if(lambda2[3] == 0 & lambda2[4] < 0){
        infimo2 <- lambda2[1] - 1/lambda2[2]
        sup2 <- Inf
      }
      infimo <- min(infimo1,infimo2)
      sup <- max(sup1,sup2)
    }
  }

  else{
    if(param == "fmkl"){
      if(lambda1[3] > 0)
        infimo <- lambda1[1] - (1/(lambda1[2] * lambda1[3]))
      if(lambda1[3] <= 0)
        infimo <- -Inf
      if(lambda1[4] > 0)
        sup <- lambda1[1] + (1/(lambda1[2] * lambda1[4]))
      if(lambda1[4] <= 0)
        sup <- Inf
    }
    if(param == "rs"){
      if(min(lambda1[3],lambda1[4]) > 0){
        infimo <- lambda1[1] - 1/lambda1[2]
        sup <- lambda1[1] + 1/lambda1[2]
      }
      if(lambda1[3] > 0 & lambda1[4] == 0){
        infimo <- lambda1[1]
        sup <- lambda1[1] + 1/lambda1[2]
      }
      if(lambda1[3] == 0 & lambda1[4] > 0){
        infimo <- lambda1[1] - 1/lambda1[2]
        sup <- lambda1[1]
      }
      if(max(lambda1[3],lambda1[4]) < 0){
        infimo <- -Inf
        sup <- Inf
      }
      if(lambda1[3] < 0 & lambda1[4] == 0){
        infimo <- -Inf
        sup <- lambda1[1] + 1/lambda1[2]
      }
      if(lambda1[3] == 0 & lambda1[4] < 0){
        infimo <- lambda1[1] - 1/lambda1[2]
        sup <- Inf
      }
    }
  }
  if(lambda0 > 0){
    sup <- max(sup,0.01)
    infimo <- min(infimo,-0.01)
  }
  return(c(infimo,sup))
}

#KS TEST
ks.test.hgld <- function(data,lambda1,lambda2,p,mixture,param,alpha,len,no.test,trace){
  if(trace)
    pb <- utils::txtProgressBar(min = 0, max = no.test, style = 3)
  lambda1 <- c(0,lambda1)
  pvalue <- vector()
  data <- data[data != 0]
  for(i in 1:no.test){
    s <- unique(sample(x = data,size = len,replace = FALSE))
    t <- ks.test(x = s,y = "phgld",lambda1 = lambda1,lambda2 = lambda2,p = p,mixture = mixture,param = param)
    pvalue[i] <- t$p.value > alpha
    if(trace)
      utils::setTxtProgressBar(pb,i)
  }
  if(trace)
    cat("\n")
  return(100 * sum(pvalue)/no.test)
}

