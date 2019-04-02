
#######################
# functions
#######################
data_importer <- function(package=NULL, ...){
#' @title bring in raw well data
#' 
#' @description bring in raw well data from github or R package
#' 
#' @details lorem ipsum
#' @param package import data from a specific package (NULL or 'wellwise')
#' @param ... ipsum
#' @import readr
#' @importFrom utils data
#' @export
#' @examples
#' runif(1)
  #require(readr)
  # read data and do elementary processing, take only single iteration of simulated data
  e = new.env()
  if(is.null(package)){
    cat("Reading in data from github\n")
    e$welldata <- read_csv("https://cirl-unc.github.io/wellwater/data/testdata.csv", col_types = cols(), ...)
  }
  if(!is.null(package)) {
    cat(paste0("Reading in data from ", package, " package\n"))
    data('welldata', package=package, envir = e)
  }
  e$welldata
  }

data_reader <- function(raw, i, 
                        expnm = c("Arsenic", "Manganese", "Lead", "Cadmium", "Copper"),
                        p=10, dx=5, N=NULL,
                        verbose=FALSE,
                        ...){
#' @title Create data in form usable for Bayesian models
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param raw data frame with raw data 
#' @param i iteration of raw data to read
#' @param expnm character vector of exposure names
#' @param p number of beta parameters in outcome model (optional - R makes a guess)
#' @param dx number of columns in data matrix (optional - R makes a guess)
#' @param N Sample size (optional - R makes a guess)
#' @param verbose print extra debugging info? FALSE
#' @param ... ipsum
#' @export
#' @importFrom stats complete.cases
#' @examples
#' runif(1)
  # read data and do elementary processingf, take only single iteration of simulated data
  if(verbose) cat(paste0("Using exposures: ", expnm, "\n"))
  #dat <- raw %>%
  #  filter(iter==i) %>%
  #  select(
  #    c("iter", "y", expnm)
  #  ) %>%
  #  filter(complete.cases(.))
  if(is.data.frame(raw)){
    dat <- raw[which(raw$iter == i), c("iter", "y", expnm), drop=FALSE]
  } else{
    raw <- as.data.frame(raw)
    dat <- raw[which(raw$iter == i), c("iter", "y", expnm), drop=FALSE]
  }
  dat <- dat[complete.cases(dat),]
  if(is.null(p)){
    p = length(expnm)
    if(verbose) cat("p is not specified, defaulting to the number of exposures\n")
  }
  if(is.null(dx)){
    dx = length(expnm)
    if(verbose) cat("dx is not specified, defaulting to the number of exposures\n")
  }
  if(is.null(N)){
    N = length(dat$y)
    if(verbose) cat("N is not specified, defaulting to the number of observed values of y\n")
  }
  
  with(dat, 
       list(y=y, 
            X=as.matrix(dat[, expnm, drop=FALSE]), 
            dx=dx, 
            p=p,
            N=N)
  )
}

get_model <- function(filename='logistic.stan', ...){
#' @title Get an existing model from the wellwise package
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param filename character string with stub, or name (stub.stan) of file in /inst directory that contains a stan model
#' @param ... UNUSED
#' @importFrom readr read_file
#' @export
#' @examples
#' runif(1)
 if(length(grep('.stan', filename))==0){
   filename=paste0(filename, ".stan")
   f = system.file('stan', filename, package='wellwise')
 }
 if(length(grep('.jags', filename))==0){
   filename=paste0(filename, ".jags")
   f = system.file('jags', filename, package='wellwise')
 }
 read_file(f)
}


convergence_check <- function(fit, ...){
#' @title Check stan model fit for convergence (not yet implemented)
#' 
#' @description Nothing here yet
#' 
#' @details Nothing here yet
#' @param fit unused
#' @param ... unused
#' @examples
#' runif(1)
  #showMethods("summary")
  #showMethods(class="stanfit")
}

data_analyst <- function(i, 
                         raw=raw, 
                         fl="~/temp.csv", 
                         s.code = NULL,
                         s.file = NULL,
                         iter=2500, 
                         warmup=500, 
                         chains=4, 
                         p = NULL,
                         dx = NULL,
                         N = NULL,
                         verbose=FALSE,
                         type='stan',
                         stanfit=NA,
                         basis = c("identity", "iqr", "sd"),
                         ...
                ){
#' @title Fit a Stan/JAGS model in a single data set
#' 
#' @description Fit a Stan/JAGS model in a single data set
#' 
#' @details Fit a Stan/JAGS model in a single data set
#' @param i iteration in simulated data
#' @param raw data frame with raw data
#' @param fl Filename to output model results for single run (useful for debugging)
#' @param s.code name of a character string with a stan model 
#' @param s.file (experimental) character name of a presepecified stan model ("wwbd_simtemplate_20190205", "wwbd_simlogistic_TransParameter_0211201920190211")
#' @param iter Number of "sweep" iterations, or iterations after the burnin/warmup
#' @param warmup warmup (stan) / n.adapt (jags). Number of iterations to allow
#' for adaptation of MCMC parameters/burnin
#' @param chains Number of parallel MCMC chains (default 4)
#' @param p number of beta parameters in outcome model (optional - R makes a guess)
#' @param dx number of columns in data matrix (optional - R makes a guess)
#' @param N Sample size (optional - R makes a guess)
#' @param verbose print extra debugging info? FALSE
#' @param type Which type of model is it? Can be 'stan' or 'jags' 
#' @param stanfit Name of stan_model or stan output that can be recycled. This
#' will use a pre-compiled version of the stan code which cuts simulation time
#' significantly over multiple runs.
#' @param basis how to treat exposure variables: 
#' "identity" = do not transform (default); 
#' "iqr"= divide every exposure variable by its interquartile range
#' "sd" = divide every exposure variable by its sd
#' Note that the parameterization of the BGF interventions will be invariant
#' to these transformations, so they change the strength of the priors, but
#' no modifications are necessary to examine interventions of the type:
#' "reduce exposure j by p%" - interventions that depend on the observed value
#' of exposure require that the basis be "identity"
#' @param ... arguments to stan() or coda.samples() for a Stan or JAGS model, 
#' respectively
#' @export
#' @import rstan rjags parallel doParallel foreach
#' @importFrom  stringr str_length
#' @importFrom  coda as.mcmc.list
#' @importFrom  readr write_csv
#' @examples
#' runif(1)
# TODO: implement garbage collection robust implementation of Stan such that
#  the model fit is preserved across iterations
  # do single analysis of data
  #stan model
  ch = NULL
  ncores = getOption("mc.cores")
  if(is.null(s.code) & is.null(s.file)) {
    if(verbose) cat("no stan model given, defaulting to logistic model with normal priors")
    s.code <- get_model('logistic')
  }
  if(str_length(s.code)<40) s.code <- get_model(s.code) # open existing model if none exists
  #make a guess at p if it's NULL and program is structured properly
  if(is.null(p) & gregexpr("beta\\[[0-9]+\\]", s.code)[[1]][1]>0) p=max(as.numeric(gsub('beta|\\[|\\]', '', # remove extraneous texts
                                       unlist(regmatches(s.code, #extract all matches
                                                         gregexpr("beta\\[[0-9]+\\]", s.code) # find all matches
                                                         ))
                                       )))
  
  sdat = data_reader(raw, i, p=p, dx=dx, N=N)
  # if basis is set to something other than 'identity' then perform a variable
  # transformation on the exposures - helpful for putting model coefficients on
  # a predictable scale
  if(basis[1] != "identity"){
    if(basis[1] == "sd"){
      sdat$X = 
        sweep(sdat$X, 2, 
            STATS=apply(sdat$X, 2, sd), FUN='/')
    }
    if(basis[1] == "range"){
      range = function(x,...) as.numeric(max(x,...) - min(x, ...))
      sdat$X = 
        sweep(sdat$X, 2, 
              STATS=apply(sdat$X, 2, range), FUN='/')
    }
    if(basis[1] == "inner90"){
      iqr = function(x,...) as.numeric(quantile(x, .95, ...) - quantile(x, .05, ...))
      if(any(apply(sdat$X, 2, iqr)==0)) 
        stop('Inner 90% change for at least one exposure variable is 0 - 
change basis to "identity", "sd", or "range"')
      sdat$X = 
        sweep(sdat$X, 2, 
              STATS=apply(sdat$X, 2, iqr), FUN='/')
    }
    if(basis[1] == "iqr"){
      iqr = function(x,...) as.numeric(quantile(x, .75, ...) - quantile(x, .25, ...))
      if(any(apply(sdat$X, 2, iqr)==0)) 
        stop('IQR for at least one exposure variable is 0 - change basis to 
"identity", "sd", "inner90", or "range"')
      sdat$X = 
        sweep(sdat$X, 2, 
              STATS=apply(sdat$X, 2, iqr), FUN='/')
    }
    if(basis[1] == "center"){
      sdat$X = 
        sweep(sdat$X, 2, 
              STATS=apply(sdat$X, 2, mean), FUN='-')
    }
    if(!(basis[1] %in% c("sd", "iqr", "inner90", "center", "range"))){
      stop('basis parameter in data_reader needs to be one of
          "sd", "iqr", "center", "range" 
          (see help file for data_analyst)')
    }
  }
  # note: keep warming up even for a while after apparent convergence: helps improve efficiency
  #  of adaptive algorithm
  if(!is.null(s.file)){
    if(type=='stan'){
     res = stan(file = system.file("stan", paste0(s.file,".stan", package = "wellwise")), 
             fit=stanfit,
             data = sdat, 
             chains = chains, 
             sample_file=fl, 
             iter=iter+warmup, 
             warmup=warmup, 
             ...)
    }
    if(type=='jags'){
        tf = system.file("jags", paste0(s.file,".jags", package = "wellwise"))
    }
  } else if(!is.null(s.code)){
    if(type=='stan'){
         res = stan(model_code = s.code,
             fit=stanfit,
             data = sdat, 
             chains = chains, 
             sample_file=fl, 
             iter=iter+warmup, 
             warmup=warmup, 
             ...)
    }
    if(type=='jags'){
      tf=tempfile()
      cat(s.code, file=tf)
    }
  }
  if(type=='jags'){
    sds = parallel.seeds("base::BaseRNG", chains);
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    res <- foreach(ch=1:chains, .packages=c('rjags')) %dopar%{
        mod = jags.model(file = tf, 
                         data = sdat, 
                         n.chains=1,
                         n.adapt=warmup,
                         inits=sds[[ch]])
        resl = coda.samples(mod,
                          variable.names='rd', # TODO: expand?
                          n.iter=iter, 
                          ...)
        resl[[1]] # outputs list by default, but just need MCMC object
    }
    stopCluster(cl)
    res = as.mcmc.list(res)
    write_csv(path=fl, as.data.frame(as.matrix(res)))
  }
  #step to include here: automated monitoring for convergence and continued sampling if not converged
  #class(res) <- 'bgfsimmod'
  # inherets 'mcmc.list' or 'stanfit' class
  res
}


sink.reset <- function(...){
#' @title clean up problems with sinking where output doesn't show up
#' 
#' @description analytic functions in wellwise use `sink` to avoid 
#' polluting the output with progress bars. A side effect is that sometimes
#' output can be turned off altogether. Runk sink.reset() after data_analyst
#' as a way to fix this problem
#' 
#' @details lorem ipsum
#' @param ... unused
#' @export
#' @examples
#' sink.reset()    
    for(i in seq_len(sink.number())){
        sink(file=NULL)
    }
}


analysis_wrapper <- function(simiters, 
                             rawdata, 
                             dir="~/temp/", 
                             root, 
                             debug=FALSE, 
                             verbose=FALSE,
                             type='stan',
                             stanfit=NA,
                             ...
                             ){
#' @title User level function to run multiple simulation iterations
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param simiters e.g. 1:2
#' @param rawdata RAW
#' @param dir output directory for intermediate files, which can usually be ignored "~/temp/"
#' @param root "test"
#' @param debug FALSE
#' @param verbose FALSE
#' @param type 'jags' or 'stan'
#' @param stanfit NA or name of object containing a fitted stan model (re-uses
#' the compiled stan model, which will be efficient for restarting iterations after
#' checking first one(s) for diagnostics
#' @param ... arguments to data analyst
#' @export
#' @examples
#' runif(1)

  call = match.call()
  # perform multiple analyses
  if(length(simiters)==1) {
    #sq = 1:simiters
    sq = simiters
  } else{
    sq = simiters
  }
  #mf = match(c("fl"), names(call))
  outfile = paste0(dir, root, "_res.csv")
  if(debug){
    outfile = "samples.csv"
    cat(paste0("Outputting samples from stan to ", outfile))
  }
  if(verbose) cat(paste0("Analyzing ", length(sq), " iterations of data\n"))
  if(verbose) cat(paste0("R output can be seen at ", paste0(dir, root, "_rmsg.txt"), "\n"))
  res = list(1:length(sq))
  j=1
  filenm = file(paste0(dir, root, "_rmsg.txt"), 'w')
  sink(filenm, split = FALSE, type = c("output", "message"))
  cat(paste0("Analyzing ", length(sq), " iterations of data\n"))
  for(i in sq){
    cat(".")
    if(type=='stan'){
      res[[j]] = data_analyst(i, rawdata, fl=outfile, verbose=verbose, type=type, 
                              stanfit=stanfit, ...)
      if(j==1) stanfit <- res[[1]] # safe file for later
    }
    # TODO: ADD type='gibbs' (alex to send example of gibbs sampler)
    if(type=='jags'){
      res[[j]] = data_analyst(i, rawdata, fl=outfile, verbose=verbose, type=type, ...)
    }
    j=j+1
  }
  #
  sink.reset()
  close(filenm)
  if(verbose) cat(paste("\nLog in ", filenm))
  class(res) <- 'bgfsimmodlist'
  res
}

summary.bgfsimmodlist <- function(
  object=NULL,
  ...
  ){
#' @name summary
#' @title Summarize results of multiple simulations
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param object ipsum
#' @param ... ipsum
#' @import rstan
#' @method summary bgfsimmodlist
#' @import rstan
#' @importFrom stats sd quantile
#' @export
#' @return a 'bgfsimres' object
#' @examples
#' runif(1)
  # summarize
  # to do
  # if have values stored from analysis_wrapper, then use those
  # else read csv files
  type = class(object[[1]])
  numf = length(object)
  postmeans = list(numf)
  divergents = numeric(numf)
  for(r in 1:numf){
    if(type=='stanfit'){
      tt = get_sampler_params(object[[r]])
      divergents[r] = sum(c(lapply(tt, function(x) sum(x[,"divergent__"]))>0))
      sum = summary(object[[r]])$summary
      idx = grep('rd', rownames(sum))
      postmeans[[r]] = sum[idx,1]
    } else{
      divergents[r] = NA
      sum = coda:::summary.mcmc.list(object[[r]])$statistics
      idx = grep('rd', rownames(sum))
      postmeans[[r]] = sum[idx,1]
    }
  }
  pm = as.data.frame(t(simplify2array(postmeans)))
  # report on posterior mean risk difference for each iteration
  ret = list(postmeans=pm,
             divergents=divergents,
    summary=data.frame(
      mean=apply(pm, 2, mean),
      sd = apply(pm, 2, sd),
      min = apply(pm, 2, min),
      max = apply(pm, 2, max),
      p25 = apply(pm, 2, function(x) as.numeric(quantile(x,p=0.25))),
      p75 = apply(pm, 2, function(x) as.numeric(quantile(x,p=0.75)))
    ))
  class(ret) <- 'bgfsimres' #defines this as an S3 object which allows us to do some shortcuts for methods
  ret
}

print.bgfsimres <- function(x, ...){
#' @name print
#' @title Print function for summarized simulation results
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param x ipsum
#' @param ... ipsum
#' @export
#' @method print bgfsimres
#' @examples
#' runif(1)
 cat('Iterations with at least one divergent chain: ', sum(x$divergents>0), '\n')
 print(x$summary, ...)
}

plot.bgfsimres <- function(x, type=ifelse(nrow(x$postmeans)>50, "density", "hist"), ...){
#' @name plot
#' @title Plot function for summarized results
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param x ipsum
#' @param type "density" "hist"
#' @param ... ipsum
#' @import ggplot2
#' @importFrom stats density
#' @importFrom tidyr gather
#' @method plot bgfsimres
#' @export
#' @examples
#' runif(1)
#' 
 #plot(density(x$postmeans), ...)
 plot_data <- gather(x$postmeans, "parm", "value", 1:ncol(x$postmeans), factor_key=TRUE)
 pbase <- ggplot(plot_data, aes_string(x = "value", color="parm")) + 
        xlab("Value") + theme_classic()
 if(length(grep("dens", type)))
   pbase + geom_density() + scale_color_discrete(name="Parameter")
 if(length(grep("hist", type)))
   pbase + geom_histogram(aes_string(fill="parm")) + 
   scale_color_discrete(name="Parameter") + 
   scale_fill_discrete(name="Parameter")
}