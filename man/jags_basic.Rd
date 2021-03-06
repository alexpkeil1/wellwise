% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{jags_basic}
\alias{jags_basic}
\title{helper: make a basic jags model}
\usage{
jags_basic(x = c("x", "z"), z = c("bmi"), y = "y", binvars = NULL,
  matx = "X", standardizex = TRUE, binary = TRUE,
  xintv = rbind(c(0.99, 0), c(0, 0.99), c(0.99, 0.99)))
}
\arguments{
\item{x}{intervenable exposures (character vector)}

\item{z}{covariates (character vector)}

\item{y}{outcome}

\item{binvars}{= non-outcome variables that are binary (character vector)}

\item{matx}{optional, name of matrix with intervenable exposures}

\item{standardizex}{logical, should x be standardized?}

\item{binary}{logical, is outcome binary?}

\item{xintv}{matrix with ncol = number of intervenable exposures, nrow = number of interventions. Each value is on [0,1] and represents the proportional decrease in the value of x upon hypothetical intervention}
}
\description{
not used by user, typically
}
\details{
lorem ipsum
}
\examples{
# library(rjags)
 dgm <- function(N=100, trueRD=0.2){
   x1 = rbinom(N, 1, 0.5)
   py00 = runif(N)*0.1 + 0.4
   l2 = rbinom(N, 1, 1/(1+exp(-1 + x1 + py00)))
   x2 = rbinom(N, 1, 1/(1+exp(-1 + x1 + l2)))
   py = py00 + trueRD*((x1 + x2)/2) #true risk difference per unit exposure;
   y = rbinom(N, 1, py)
   data.frame(x1, l2, x2, y)
 }
 dat = as.list(dgm(100))
 dat$N = 100
 dat$p = 5
 
 source("~/Epiprojects/wellwater/sims/code/make_stan_terms.R")
 
 mod = jags_basic(x=c('x1', 'x2'), z = 'l2', y='y', 
   binvars=c('x1', 'x2', 'l2'), xintv = rbind(c(1,0), 
   c(0,1),c(1,1)), binary=TRUE, matx = NULL)
 cat(mod)
# usage in jags (or edit by hand)
# not run
# tf = tempfile()
# cat(mod, file=tf)
# jags.model(file = tf, data = dat, n.chains=1)
}
