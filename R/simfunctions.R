
#######################
# functions
#######################
data_importer <- function(package=NULL, fast=TRUE,...){
#' @title bring in raw well data
#' 
#' @description bring in raw well data from github or R package
#' 
#' @details lorem ipsum
#' @param package import data from a specific package (NULL or 'wellwise')
#' @param fast logical, if TRUE (default), downloads a smaller data set for faster testing. 
#' If FALSE, downloads more realistic data for use in simulations.
#' @param ... ipsum
#' @importFrom utils data read.csv download.file
#' @export
#' @examples
#' runif(1)
# #' @import readr
  #require(readr)
  # read data and do elementary processing, take only single iteration of simulated data
  e = new.env()
  if(is.null(package)){
    cat("Reading in data from github\n")
    cat("try data_importer(package='wellwise') to speed this up substantially!\n")
    #e$welldata <- read_csv("https://cirl-unc.github.io/wellwater/data/testdata20190404.csv", col_types = cols(.default = col_double()), ...)
    if(fast) e$welldata <- read.csv("https://cirl-unc.github.io/wellwater/data/testdata.csv")
    if(!fast) {
      t = tempfile()
      download.file("https://cirl-unc.github.io/wellwater/data/welldata.rda", t)
      e$welldata <- load(t)
    }
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


julia.model <- function(j.code,
                       iter,
                       warmup, 
                       chains, 
                       sdat=sdat,
                       fl=tempfile(),
                       jinstance=NULL,
                       cleanafter=FALSE,
                       verbose=TRUE,
                       debug=FALSE,
                       addpackages=TRUE,
                       juliabin=NULL
                       ){
#' @title Fit a julia gibbs sampler (given by s.code)
#' 
#' @description Fit a julia model to a single data set
#' 
#' @details Fit a julia model to a single data set
#' @param j.code name of a character string with a stan model 
#' @param iter Number of "sweep" iterations, or iterations after the burnin/warmup
#' @param warmup warmup (stan) / n.adapt (jags) / burnin (julia). Number of iterations to allow
#' for adaptation of MCMC parameters/burnin
#' @param chains Number of parallel MCMC chains (default 4)
#' @param sdat list of data created from 
#' @param fl Filename to output model results for single run (useful for debugging)
#' @param jinstance existing julia_setup() instance
#' @param cleanafter should the extra workers be killed after running?
#' @param verbose print extra stuff 
#' @param debug print extra debugging info? FALSE
#' @param addpackages install/reinstall julia packages? TRUE (default) is slower, but should be
#' done the first time running each model. Afterwards, can be set to FALSE to save
#' a good chunk of computing time.
#' @param juliabin NULL or path to directory containing the julia binary file
#'  e.g. /Applications/Julia-1.1.app/Contents/Resources/julia/bin/ on mac
#'  C:/Users/<username>]/AppData/Local/Julia-1.1.0/bin/ on windows 10
#' @import JuliaCall
#' @importFrom  coda as.mcmc.list mcmc
#' @export
  call = match.call()
  print(call)
  if(verbose) cat(paste0("Julia code stored in ", fl))
  runonce <- c(
    "print(\"Testing wellwisejl\")"
  )
  runworker <- c(
     "using Pkg, Distributed, LinearAlgebra, Statistics, Distributions,
      DataFrames, PolyaGammaDistribution, AltDistributions, wellwisejl",
     j.code,
    'println("julia utility commands done")'
  )
  cat("# Julia worker command file created by wellwise package \n\n\n", file=fl)
  for(u in runworker){
    if(debug) cat(u, "\n\n")
    cat(u, "\n\n", file=fl, append=TRUE)
  }

  #cleanup commands
  cleanup <- c(
    "for p in procs() 
      if p>1 
       rmprocs(p)
      end
    end",
    "println(\"cleaning commands finished\")"
  )
  # prelims
  if(is.null(jinstance)) {
    if(is.null(juliabin)) julia <- julia_setup()
    if(!is.null(juliabin)) julia <- julia_setup(JULIA_HOME = juliabin)
    julia$command("include(x) = Base.include(Main, x)") # hack because include doesnt work with embedded julia
    pkgs <- c("Pkg", "Distributed", "LinearAlgebra", "Statistics", "Distributions",
             "DataFrames")
    pkgs2 <- c("https://github.com/alexpkeil1/wellwisejl",
               "https://github.com/alexpkeil1/PolyaGammaDistribution.jl", 
             "https://github.com/alexpkeil1/AltDistributions.jl"
             )
    if(addpackages){
      for(pkg in pkgs) {
        print(pkg)
        julia$install_package_if_needed(pkg)
      }
      for(pkg in pkgs2) {
        print(pkg)
        julia$install_package(pkg)
      }
    }
    # create some necessary functions  
    cat("\n", file=normalizePath(paste0(fl, '2'), mustWork = FALSE), append=FALSE)
    for(u in runonce) {
      if(debug) print(u)
      cat(u, "\n\n", file=normalizePath(paste0(fl, '2')), append=TRUE)
    }
    #  call the "runonce" code
    julia_source(normalizePath(paste0(fl, '2'), winslash = "/"))
    #julia$command(paste0("include(\"",normalizePath(fl, winslash = "/"),"2\")"))
    # export commands to workers, this requires slightly different path naming on windows
    #  due to order of evaluation vs. escaping
    julia$command(paste0("include(\"",normalizePath(fl, winslash = "/"),"\")"))
    julia$command(paste0("addmoreprocs(", chains, ")"))
    julia$command(paste0("@everywhere include(\"",normalizePath(fl, winslash = "/"),"\")"))
    } else {
      if(verbose) cat("Utilizing existing instance of julia_setup()")
      julia <- jinstance
    }
  # read in data
  julia$assign("rdat", data.frame(cbind(y=sdat$y, sdat$X)))

  if(verbose) cat(paste(readLines(normalizePath(fl)), collapse = "\n"))
  
  # run model on all workers and collect results
  jcmd = paste0('r = runmod(gibbs, rdat, ',iter+warmup,', ',warmup, ', ', chains,')')
  if(verbose) print(paste("running models with '", jcmd, "'"))
  julia$command(jcmd) # run gibbs sampler in Julia
  jres = julia$eval('r') # bring julia results into R as matrix
  res = as.mcmc.list(lapply(jres, function(x) mcmc(x, start=warmup+1, end=iter+warmup)))
  if(cleanafter){
    for(c in cleanup) {
      if(debug) print(c)
      julia$command(c) # also number of chains -todo move this to higher level control
    }
  }
  list(res=res, jinstance=julia)
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
                         env = NULL,
                         verbose=FALSE,
                         debug=FALSE,
                         type=c('stan', 'jags', 'julia'),
                         stanfit=NA,
                         basis = c("identity", "iqr", "sd"),
                         juliabin=NULL,
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
#' @param env environment for existing julia instance
#' @param verbose print extra debugging info? FALSE
#' @param debug print extra debugging info? FALSE
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
#' @param juliabin NULL or path to directory containing the julia binary file
#'  e.g. /Applications/Julia-1.1.app/Contents/Resources/julia/bin/ on mac
#'  C:/Users/<username>]/AppData/Local/Julia-1.1.0/bin/ on windows 10
#' @param ... arguments to stan(), coda.samples(), or julia.sample() for a Stan, JAGS, or julia model, 
#' respectively
#' @export
#' @import rstan rjags parallel doParallel foreach JuliaCall
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
  type = type[1]
  ncores = getOption("mc.cores")
  if(is.null(ncores)){
    print(paste0('mc.cores option not set, setting to ', chains))
    options(mc.cores = chains)
  }
  if(!is.null(ncores) && ncores< chains){
    print(paste0('mc.cores option less than number of chains, setting to ', chains))
    print(paste0(''))
    options(mc.cores = chains)
  }
  if(getOption("mc.cores") > parallel::detectCores()){
    stop(
      paste0("wellwise: not enough computer cores to run this simulation in parallel - set chains to ", 
             parallel::detectCores(), 
             " or fewer"))
  }
  if(is.null(s.code) & is.null(s.file)) {
    if(verbose) stop("no stan/jags/julia model given, (specify either s.code or s.file[experimental])")
    #s.code <- get_model('logistic')
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
    }else if(type=='jags'){
        tf = system.file("jags", paste0(s.file,".jags", package = "wellwise"))
    } else{
      stop(paste("s.file not implemented for type", type, ". Try setting s.file=NULL"))
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
    } else if(type=='jags'){
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
  if(type=='julia'){
    if(verbose) cat("julia models are experimental")
    if(fl=="") stop("fl must be a writable filename for type='julia'")
    if(!exists("jinstance", envir=env)) env$jinstance=NULL
    resj = julia.model(iter=iter,
                   j.code=s.code,
                   warmup=warmup, 
                   chains=chains, 
                   sdat=sdat,
                   fl=fl,
                   verbose=verbose,
                   jinstance = env$jinstance,
                   cleanafter = FALSE, 
                   debug=debug,
                   juliabin=juliabin,
                   ...
                  )
    res = resj[['res']]
    env$jinstance = resj[['jinstance']]
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
                             root='', 
                             debug=FALSE, 
                             verbose=FALSE,
                             type='stan',
                             stanfit=NA,
                             juliabin=NULL,
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
#' @param type 'jags', 'stan', or 'julia'
#' @param stanfit NA or name of object containing a fitted stan model (re-uses
#' the compiled stan model, which will be efficient for restarting iterations after
#' checking first one(s) for diagnostics
#' @param juliabin NULL or path to directory containing the julia binary file
#'  e.g. /Applications/Julia-1.1.app/Contents/Resources/julia/bin/ on mac
#'  C:/Users/<username>]/AppData/Local/Julia-1.1.0/bin/ on windows 10
#' @param ... arguments to data analyst
#' @importFrom utils write.csv
#' @export
#' @examples
#' runif(1)
  t1 <- Sys.time()
  call = match.call()
  # perform multiple analyses
  if(length(simiters)==1) {
    #sq = 1:simiters
    sq = simiters
  } else{
    sq = simiters
  }
  #mf = match(c("fl"), names(call))
  dir = normalizePath(path.expand(dir), mustWork=TRUE)
  outfile = normalizePath(file.path(dir, paste0(root, "_res.csv")), mustWork=FALSE)
  if(debug & type=='stan'){
    outfile = "samples.csv"
    cat(paste0("Outputting samples from stan to ", outfile))
  }
  if(verbose) cat(paste0("Analyzing ", length(sq), " iterations of data\n"))
  if(verbose) cat(paste0("R output can be seen at ", paste0(dir, root, "_rmsg.txt"), "\n"))
  res = list(1:length(sq))
  j=1
  filenm = paste0(dir, root, "_rmsg.txt")
  if(!debug) sink(filenm, split = FALSE, type = c("output", "message"))
  if(debug) sink(filenm, split = TRUE, type = c("output", "message"))
  cat(paste0("Analyzing ", length(sq), " iterations of data\n"))
  if(type == 'julia'){
    checkjulia()
  }
  for(i in sq){
    cat(".")
    if(type=='stan'){
      res[[j]] = data_analyst(i, rawdata, fl=outfile, verbose=verbose, debug=debug, 
                              type=type, stanfit=stanfit, ...)
      if(j==1) stanfit <- res[[1]] # save file for later
    }
    # TODO: ADD type='gibbs' (alex to send example of gibbs sampler)
    if(type == 'jags'){
      res[[j]] = data_analyst(i, rawdata, fl=outfile, verbose=verbose, debug=debug, 
                              type=type, ...)
    }
    if(type == 'julia'){
      jfn = normalizePath(file.path(dir, paste0(root, "_jcode.txt")), mustWork=FALSE)
      if(verbose) cat(paste0("Julia code can be seen at ", jfn, "\n"))
      jenv = new.env()
      res[[j]] = data_analyst(i, rawdata, fl=jfn, verbose=verbose, debug=debug, 
                              type=type, env=jenv, juliabin=juliabin, ...)
      write.csv(do.call("rbind", res[[j]]), outfile)
    }
    j=j+1
  }
  #
  sink.reset()
  t2 <- Sys.time()
  print(paste("Run time :", t2-t1, "mins"))
  #close(filenm)
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
      sum = summary.mcmc.list(object[[r]])$statistics
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



######### 
# helpers
##########

checkjulia <- function(juliabin=NULL){
  curr = Sys.info()["sysname"]
  joiner = ":"
  jpath = ""
  if(!is.null(juliabin)){
    jpath = c(jpath, normalizePath(juliabin))
  }
  if(Sys.which("julia") != ""){
    np = normalizePath(dirname(Sys.which("julia")))
    jpath = c(jpath, np)
  }
  jhome = Sys.getenv("JULIA_HOME")
  path = Sys.getenv("PATH")
  jpath = jpath[jpath!=""]
  if (tolower(substr(curr, 1, 3)) == "win") {
    joiner = ";"
    un = Sys.info()["user"]
    paths = c("/AppData/Local/Julia-0.7.0/bin/",
              "/AppData/Local/Julia-1.0.0/bin/",
              "/AppData/Local/Julia-1.0.1/bin/",
              "/AppData/Local/Julia-1.0.2/bin/",
              "/AppData/Local/Julia-1.0.3/bin/",
              "/AppData/Local/Julia-1.1.0/bin/")
    for(pth in rev(paths)){
      jpath = suppressWarnings(c(jpath, normalizePath(file.path("C:/Users/", un, pth))))
    }
  }
  if (Sys.getenv("JULIA_HOME") != "") {
    Sys.setenv(JULIA_HOME = paste0(c(jpath, jhome), collapse = joiner))
  }else Sys.setenv(JULIA_HOME = paste0(jpath, collapse = joiner))
  if (Sys.getenv("PATH") != "") {
    Sys.setenv(PATH = paste0(c(jpath, path), collapse = joiner))
  }else Sys.setenv(PATH = paste0(jpath, collapse = joiner))
  for(j in jpath){
    options("JULIA_HOME"=j)
    jp = suppressWarnings(JuliaCall:::julia_locate())
    if(!is.null(jp)) break
  }
  
  res = try(JuliaCall::julia_setup(JULIA_HOME = jp, useRCall = FALSE, 
                                   force = TRUE, verbose = FALSE), silent = TRUE)
  if (class(res) == "try-error") {
    ret = "Julia is not working, try setting 'juliabin' parameter to path containing julia: e.g. C:/Users/<username>/AppData/Local/Julia-1.1.0/bin/ on windows 10"
  }else ret = "Julia appears to be working"
  ret
}



summary.mcmc.list <- function (...) 
{
#' @name summary.mcmc.list
#' @title Summary method for mcmc.list object from coda package
#' 
#' @description This function was taken directly from the `coda` package. Not
#' exported in that package.
#' 
#' @details lorem ipsum
#' @param ... arguments to coda:::summary.mcmc.list (object, quantiles=c(.025, .25, .5, .75, .975))
#' @importFrom utils getFromNamespace
#' @method summary mcmc.list
#' @export
#' @examples
#' runif(1)
  fun = getFromNamespace("summary.mcmc.list", "coda")
  fun(...)
}
