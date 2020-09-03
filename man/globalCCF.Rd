\name{globalCCF}
\alias{globalCCF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Global Cross Covariance/Correlation Function Estimation for Multiple Markov Chains
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates the global cross-correlation (or covariance) function (G-CCF) of two component for a particular chain from a given list of parallel Markov chains sampled using same Markov transition kernel. The function differs from the base `ccf` function by allowing the users to center the chains around the global mean from all Markov chains. A different mean argument can also be specified by the user.
}
\usage{
globalCCF(x, chain = NULL, components = NULL, lag.max = NULL, type = "correlation", 
  mean = "global", plot = TRUE, ymin = NULL, ymax = NULL, main = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    a list of numeric \eqn{n \times p} matrices representing parallel Markov chains where \eqn{n} is sample size and \eqn{p} is dimension.
  }
  \item{chain}{
    numeric. The chain for which CCF is to be calculated.
  }
  \item{components}{
  a numeric vector of two components between which CCF is to be calculated.
  }
  \item{lag.max}{
  maximum lag for ACF calculations. Default is determined by base \code{ccf}. 
  }
  \item{type}{
  character string giving the type of CCF to be computed. Input ``covariance'' for calculating cross-covariance function and ``correlation'' for calculating cross-correlation function.
  }
  \item{mean}{
    character string or vector of length \eqn{p}. Tells how to center the chain for CCF calculations. Allowed character strings are ``local'' or ``global''. User can specify a mean vector of length \eqn{p}.
  }
  \item{plot}{
    logical. If \code{TRUE}, plot the CCF returned by the function.
  }

  \item{ymin}{
    numeric. If \code{plot == TRUE}, sets the lower limit of y scale.
  }
  \item{ymax}{
    numeric. If \code{plot == TRUE}, sets the upper limit of y scale.
  }
  \item{main}{
  character string. If \code{plot == TRUE}, sets the main title on top of the plot.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The function calculates sample CCF between the given two components; but for only one chain from a list of parallel Markov chains. \code{mean} parameter centers the concerned chain around its sample mean ("local"), overall sample mean of all chains ("global") or any user specified mean (say, when the true mean is known).

}
\value{
An object of class ``acf''.
}

%\references{
%% ~put references to the literature/web site here ~
%}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{globalACF}} which gives ACF for multiple Markov chains.
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

N <- 5e2 # run for smaller length due to time
p <- 0.7
mu1 <- c(-4,-4)
mu2 <- c(4,4)
sigma1 <- diag(2)
sigma2 <- diag(c(1.5,1.5))

chain1 <- mh.mcmc(start = c(-4,-4), p, mu1, mu2, sigma1, sigma2, N, 1)
chain2 <- mh.mcmc(start = c(4,4), p, mu1, mu2, sigma1, sigma2, N, 1)
x <- list(as.matrix(chain1), as.matrix(chain2)) # Markov chain list

par(mfrow = c(1,2))
# bad CCF. Same as ccf()
globalCCF(x, chain=1, type="correlation", mean="local", ymin=0)  # bad CCF
globalCCF(x, chain=1, type="correlation", mean="global", ymin=0) # good CCF
}

\references{

  Agarwal, M and Vats, D (2020+)  Globally-centered autocovariances in MCMC, \emph{arxiv}.
}
