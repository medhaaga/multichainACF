
globalACF <- function(x, chains = NULL, component = NULL, lag.max = NULL,
                      mean = "global", type = "correlation",
                      avg = TRUE, plot = TRUE, leg = FALSE, col = "black", ymin = NULL, ymax = NULL, main = "")
{
  if(is.null(chains) && avg == FALSE)
    stop("Nothing to return/plot")

  if(class(x) != "list")
    stop("a list of Markov chains is required")

  for (i in 1:length(x)){
    if(!is.numeric(x[[i]]))
      stop("'x' must be numeric")
  }

  m <- length(x)
  n <- as.integer(nrow(x[[1]]))
  p <- as.integer(ncol(x[[1]]))

  if (is.na(n) || is.na(p))
    stop("no. of samples and dim must be integer")

  if(p > 1 && is.null(component))
    stop("component parameter missing")

  if(component > p)
    stop("Incorrect component")

  if(p==1 && is.null(component))
    component = 1

  if (is.null(lag.max))
    lag.max <- floor(10 * (log10(n)))
  lag.max <- as.integer(min(lag.max, n - 1L))

  if (is.na(lag.max) || lag.max < 0)
    stop("'lag.max' must be at least 0")

  if (!is.null(chains) && length(chains) == 1 && chains == 0)
    chains <- seq(1,m,1)

  if(!is.null(chains))
    l <- length(chains)

  if((!is.null(ymin) && !is.numeric(ymin)) || (!is.null(ymax) && !is.numeric(ymax))){
    warnings("ylim is not correct. Default ylim will be used.")
    ymin=NULL
    ymax=NULL
  }

  mc.chains <- array(0, dim = c(n, m))
  for (i in 1:m)
    mc.chains[,i] <- x[[i]][, component]

  ################## centering ######################
  
  if(!is.numeric(mean) && mean == "global"){
    center = mean(mc.chains)
    mc.chains <- mc.chains - center
  } else if (!is.numeric(mean) && mean == "local"){
    mc.chains <- scale(mc.chains, scale = FALSE)
  } else {
    if(!is.numeric(mean) || length(mean) != p)
      stop("'mean' is not numeric or of correct dimensions.")
    mc.chains <- apply(mc.chains, 2, '-', mean[component])
  }
  #######################################################

  acf.list <- list()
  avg.acf <- list()
  attr(avg.acf, "class") <- "acf"

  if (avg)
  {
    acf.list[[1]] <- acf(mc.chains[,1], type = type, plot = FALSE, demean = FALSE, lag.max = lag.max)
    avg.acf <- acf.list[[1]]

    for (i in 2:m){
      acf.list[[i]] <- acf(mc.chains[,i], type = type, plot = FALSE, demean = FALSE, lag.max = lag.max)
      avg.acf$acf <- avg.acf$acf + acf.list[[i]]$acf
    }

    avg.acf$acf <- avg.acf$acf/m

    if (plot)
    {
      if(is.null(ymin))
        ymin <- min(unlist(lapply(acf.list, function(x) min(x$acf))))

      if(type == "correlation")
        ymin <- min(0, ymin)

      if (is.null(ymax))
        ymax <- max(unlist(lapply(acf.list, function(x) max(x$acf))))

      ylim <- c(ymin, ymax)


      if(!is.null(chains) && length(chains)>1)
        {
          plot(avg.acf, ci=0, type = "l", col = col,
               lwd = 2, xlab = "lag", ylab = "ACF", ylim = ylim, main = main)

          for (i in 1:l){
            j <- chains[i]
            lines(0:lag.max, as.matrix(acf.list[[j]]$acf), type = "l", col = adjustcolor(col, alpha.f = .5), lwd = 2, lty = 2, yaxt = 'n', xaxt = 'n')
          }

          if(leg)
            legend("bottomleft", legend = c("Avg-ACF", "ACF"),
                 col = c(col, adjustcolor(col, alpha.f = 0.5)),
                 lty = c(1,2), lwd = 2)
        }  else{

          if(!is.null(chains) && length(chains)==1){
            plot(acf.list[[chains]], type = 'h', xlab = "Lag", ylab = "ACF", ylim = ylim, main = main)
            lines(0:lag.max, avg.acf$acf, type = "l", col = col, lwd = 2, lty = 1, yaxt = 'n', xaxt = 'n')
            if(leg)
              legend("bottomleft", legend = "Avg-ACF", col = col, lty = 1, lwd = 2)
          } else
            plot(avg.acf, ci=0, xlab = "Lag", ylab = "ACF", ylim = ylim, main = main)
        }
    }
  } else
  {

    for (j in 1:l){
      acf.list[[j]] <- acf(mc.chains[,chains[j]], type = type, plot = FALSE, demean = FALSE, lag.max = lag.max)
    }

    if(plot)
    {
      if(is.null(ymin))
        ymin <- min(unlist(lapply(acf.list, function(x) min(x$acf))))

      if(type == "correlation")
        ymin <- min(0, ymin)

      if (is.null(ymax))
        ymax <- max(unlist(lapply(acf.list, function(x) max(x$acf))))

      ylim = c(ymin, ymax)

      if (length(chains) == 1)
        plot(acf.list[[1]], xlab = "Lag", ylab = "ACF", ylim = ylim, main = main) else
        {
          plot(as.matrix(acf.list[[1]]$acf), type = "l", col = adjustcolor(col, alpha.f = .5), lwd = 2, lty=1, xlab = "lag", ylab = "ACF", ylim = ylim, main = main)
          for (i in 2:l)
          {
            lines(0:lag.max, as.matrix(acf.list[[i]]$acf), type = "l", col = adjustcolor(col, alpha.f = .5), lwd = 2, lty=1)
          }
        }
    }
  }
  if(avg)
    invisible (list("avgACF" = avg.acf, "chainsACF" = acf.list)) else
      invisible(acf.list)
}



globalCCF <- function(x, chain = NULL, components = NULL, lag.max = NULL,
                      type = "correlation", mean = "global",
                      plot = TRUE, ymin = NULL, ymax = NULL, main = ""){

  if(is.null(chain))
    stop("'chain' parameter missing")

  if(class(x) != "list")
    stop("a list of Markov chains is required")

  for (i in 1:length(x)){
    if(!is.numeric(x[[i]]))
      stop("'x' must be numeric")
  }

  n <- as.integer(nrow(x[[1]]))
  p <- as.integer(ncol(x[[1]]))

  if (is.na(n) || is.na(p))
    stop("no. of samples and dim must be integer")

  if (p==2 && is.null(components))
    components <- c(1,2)

  if(p <2 || is.null(components) || length(components) != 2)
    stop("two components needed")

  if (is.null(lag.max))
    lag.max <- floor(10 * (log10(n)))
  lag.max <- as.integer(min(lag.max, n - 1L))

  if (is.na(lag.max) || lag.max < 0)
    stop("'lag.max' must be at least 0")


  if((!is.null(ymin) && !is.numeric(ymin)) || (!is.null(ymax) && !is.numeric(ymax))){
    warnings("ylim is not correct. Default ylim is be used.")
    ymin=NULL
    ymax=NULL
  }

  mc.chain <- x[[chain]][,components]
  mu <- colMeans(mc.chain)

  ################## setting the center parameter ######################
  if(is.null(mean) || mean == "global"){
    center = colMeans(Reduce("rbind", x))[components]
  } else if (mean == "local"){
    center <- mu
  } else {
    if(!is.numeric(mean) || length(mean) != p)
      stop("'mean' is not numeric or of correct dimensions.")
    center <- mean[components]
  }
  #######################################################


  global.ccf <- ccf(mc.chain[,1], mc.chain[,2], type = "covariance", plot = FALSE, lag.max = lag.max)
  x.cen <- scale(mc.chain[,1], center = TRUE, scale = FALSE)
  y.cen <- scale(mc.chain[,2], center = TRUE, scale = FALSE)
  x.glob.cen <- scale(mc.chain[,1], center = center[1], scale = FALSE)
  y.glob.cen <- scale(mc.chain[,2], center = center[2], scale = FALSE)
  for (k in 1:lag.max){
    res1 <- (n - k)*(mu[1] - center[1])*(mu[2] - center[2]) - (mu[1] - center[1])*sum(y.cen[1:k]) - (mu[2] - center[2])*sum(x.cen[(n-k+1):n])
    res2 <- (n - k)*(mu[1] - center[1])*(mu[2] - center[2]) - (mu[1] - center[1])*sum(y.cen[(n-k+1):n]) - (mu[2] - center[2])*sum(x.cen[1:k])
    global.ccf$acf[lag.max+1+k] = global.ccf$acf[lag.max+1+k] + res1/n
    global.ccf$acf[lag.max+1-k] = global.ccf$acf[lag.max+1-k] + res2/n
  }
  global.ccf$acf[lag.max+1] = global.ccf$acf[lag.max+1] + (mu[1] - center[1])*(mu[2] - center[2])
  if (type == "correlation")
  {var.x <- sum(x.glob.cen * x.glob.cen)/n
  var.y <- sum(y.glob.cen * y.glob.cen)/n
  global.ccf$acf = global.ccf$acf/sqrt(var.x*var.y)}

  if (plot){

    if(is.null(ymin))
      ymin <- min(global.ccf$acf)
    
    if(type == "correlation")
      ymin <- min(0, ymin)

    if (is.null(ymax))
      ymax <- max(global.ccf$acf)

    plot(global.ccf, xlab = "Lag", ylab = "CCF", ylim = c(ymin, ymax), main = main)
    }
  invisible (global.ccf)
  }
