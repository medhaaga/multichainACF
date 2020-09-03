\name{multichainACF-package}
\alias{multichainACF}
\docType{package}
\title{
Autocorrelation function (ACF) plots for Multiple Markov Chains
}
\description{
Given a list of multiple Markov chains, this package offers functions for ACF and CCF plots  that pools information from all chains yielding more accurate estimates.
}

\details{
\tabular{ll}{
Package: \tab multichainACF\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2020-08-26\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Medha Agarwal <medhaaga@gmail.com> \\
Dootika Vats <dootika@iitk.ac.in> \\

Maintainer: Medha Agarwal <medhaaga@gmail.com>
}

\references{

  Agarwal, M and Vats, D (2020+)  Globally-centered autocovariances in MCMC, \emph{arxiv}.
}


\examples{
## Gaussian mixture target using random walk Metropolis Hastings

set.seed(1)
library(mvtnorm)
log.density <- function(x, p, mu1, mu2, sigma1, sigma2){
  return(log(p*dmvnorm(x, mean = mu1, sigma = sigma1)
    + (1-p)*dmvnorm(x, mean = mu2, sigma= sigma2)))
}


mh.mcmc <- function(start, p, mu1, mu2, sigma1, sigma2, N, h){
  X <- matrix(0, nrow = N, ncol = 2)
  X[1,] <- start
  acc <- 0
  for (i in 2:N){
    prop <- rnorm(2, mean = X[i-1,], sd = h)
    ratio <- log.density(prop, p, mu1, mu2, sigma1, sigma2)
    - log.density(X[i-1,], p, mu1, mu2, sigma1, sigma2)
    if(runif(1) < exp(ratio)){
      X[i,] <- prop
      acc <- acc + 1
    } else
      X[i,] <- X[i-1]
  }
  print (acc/N)
  return (X)
}

N <- 5e2  # run for smaller length due to time
p <- 0.7
mu1 <- c(-4,-4)
mu2 <- c(4,4)
sigma1 <- diag(2)
sigma2 <- diag(c(1.5,1.5))

chain1 <- mh.mcmc(start = c(-4,-4), p, mu1, mu2, sigma1, sigma2, N, 1)
chain2 <- mh.mcmc(start = c(4,4), p, mu1, mu2, sigma1, sigma2, N, 1)
x <- list(as.matrix(chain1), as.matrix(chain2)) # Markov chain list


par(mfrow = c(2,2))
#bad ACF same as acf()

globalACF(x, chains=0, component=1, type="correlation", mean="local") #bad ACF
globalACF(x, chains=0, component=1, type="correlation", mean="global") # good ACF

globalACF(x, chains=1, component=1, type="correlation", mean="local", avg = TRUE) #bad ACF
globalACF(x, chains=1, component=2, type="correlation", mean="global", avg = TRUE) # good ACF

}
