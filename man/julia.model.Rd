% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simfunctions.R
\name{julia.model}
\alias{julia.model}
\title{Fit a julia gibbs sampler (given by s.code)}
\usage{
julia.model(j.code, iter, warmup, chains, sdat = sdat, fl = tempfile(),
  jinstance = NULL, cleanafter = FALSE, verbose = TRUE,
  debug = FALSE, addpackages = TRUE, juliabin = NULL)
}
\arguments{
\item{j.code}{name of a character string with a stan model}

\item{iter}{Number of "sweep" iterations, or iterations after the burnin/warmup}

\item{warmup}{warmup (stan) / n.adapt (jags) / burnin (julia). Number of iterations to allow
for adaptation of MCMC parameters/burnin}

\item{chains}{Number of parallel MCMC chains (default 4)}

\item{sdat}{list of data created from}

\item{fl}{Filename to output model results for single run (useful for debugging)}

\item{jinstance}{existing julia_setup() instance}

\item{cleanafter}{should the extra workers be killed after running?}

\item{verbose}{print extra stuff}

\item{debug}{print extra debugging info? FALSE}

\item{addpackages}{install/reinstall julia packages? TRUE (default) is slower, but should be
done the first time running each model. Afterwards, can be set to FALSE to save
a good chunk of computing time.}

\item{juliabin}{NULL or path to directory containing the julia binary file
e.g. /Applications/Julia-1.1.app/Contents/Resources/julia/bin/ on mac
C:/Users/<username>]/AppData/Local/Julia-1.1.0/bin/ on windows 10}
}
\description{
Fit a julia model to a single data set
}
\details{
Fit a julia model to a single data set
}
