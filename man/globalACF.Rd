\name{globalACF}
\alias{globalACF}

\title{

Global Auto Covariance/Correlation Function Estimation for Multiple Markov Chains
}
\description{

Calculates the global autocorrelation (or covariance) function (G-ACF) of a component for a given list of multiple Markov chains sampled using the same Markov transition kernel. The function differs from the base `acf` function by allowing the users to center the chains around the global mean from all Markov chains. A different mean argument can also be specified by the user.
}
\usage{
globalACF(x, chains = NULL, component = NULL, lag.max = NULL, 
mean = "global", type = "correlation", avg = TRUE, plot = TRUE, 
leg = FALSE, col = "black", ymin = NULL, ymax = NULL, main = "")
}

\arguments{
  \item{x}{
  a list of numeric \eqn{n \times p} matrices representing parallel Markov chains where \eqn{n} is sample size and \eqn{p} is dimension of state space.
  }
  \item{chains}{
    numeric vector. The chains for which G-ACF are to be calculated and returned. Input \code{chains = 0} to return for all chains in x.
  }
  \item{component}{
  numeric. Component for which ACF is to be calculated for all the given chains.
  }
  \item{lag.max}{
  maximum lag for ACF calculations. Default is determined by base \code{acf}. 
  }
  \item{mean}{
  character string or vector of length \eqn{p}. Tells how to center all the chains for ACF calculations. Allowed character strings are ``local'' or ``global''. User can specify a mean vector of length \eqn{p}.
  }
  \item{type}{
  character string giving the type of ACF to be computed. Input ``covariance'' for calculating autocovariance function and ``correlation'' for calculating ``autocorrelation'' function.
  }
  \item{avg}{
  logical. If \code{TRUE}, the average ACF over all chains is also returned.
  }
  \item{plot}{
  logical. If \code{TRUE}, plot all the ACFs that are returned by the function.
  }
  \item{leg}{
  logical. If \code{TRUE}, show legend.
  }
  \item{col}{
  Determines the color of line ACF plots.
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
The function calculates and plots the sample ACF for any number of chains out of the given list of multiple Markov chains; but for only one component. The ``local'' and ``global'' mean centers the chains around their individual mean or the overall mean of the chains. User can also specify a mean (say, when the true mean is known).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{avgACF }{if avg is \code{TRUE}, ACF averaged over all the chains}
  \item{chainsACF }{list of ACF for all the chains provided in \code{chains} parameter}
}
%\references{
%% ~put references to the literature/web site here ~
%}

\seealso{
\code{\link{globalCCF}} which gives CCF for multiple Markov chains.
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
globalACF(x, chains=0, component=1, type="correlation", mean="local")
globalACF(x, chains=0, component=1, type="correlation", mean="global") # good ACF

globalACF(x, chains=1, component=1, type="correlation", mean="local", avg = TRUE) #bad ACF
globalACF(x, chains=1, component=2, type="correlation", mean="global", avg = TRUE) # good ACF

}

\references{

  Agarwal, M and Vats, D (2020+)  Globally-centered autocovariances in MCMC, \emph{arxiv}.
}
