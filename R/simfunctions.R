# functions
#######################
data_importer <- function(package=NULL, ...){
#' @title lorem ipsum
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param package import data from a specific package (NULL or 'wellwise')
#' @param ... ipsum
#' @importFrom readr read_csv cols
#' @export
#' @examples
#' runif(1)
  #require(readr)
  # read data and do elementary processing, take only single iteration of simulated data
  cat("Reading in data from github\n")
  welldata <- read_csv("https://cirl-unc.github.io/wellwater/data/testdata.csv", col_types = cols(), ...)
  if(!is.null(package))load("https://cirl-unc.github.io/wellwise/tree/master/data/welldata.RData")
  welldata
  }

data_reader <- function(raw, i, 
                        expnm = c("Arsenic", "Manganese", "Lead", "Cadmium", "Copper"),
                        p=10, dx=5, N=NULL,
                        ...){
#' @title lorem ipsum
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param raw ipsum
#' @param i ipsum
#' @param expnm ipsum
#' @param p ipsum
#' @param dx ipsum
#' @param N ipsum
#' @param ... ipsum
#' @export
#' @import dplyr
#' @importFrom stats complete.cases
#' @examples
#' runif(1)
  #require(dplyr)
  # read data and do elementary processingf, take only single iteration of simulated data
  cat(paste0("Using exposures: ", expnm, "\n"))
  dat <- raw %>%
    filter(iter==i) %>%
    select(
      c("iter", "y", expnm)
    ) %>%
    filter(complete.cases(.))
  if(is.null(p)){
    p = length(expnm)
    cat("p is not specified, defaulting to the number of exposures\n")
  }
  if(is.null(dx)){
    dx = length(expnm)
    cat("dx is not specified, defaulting to the number of exposures\n")
  }
  if(is.null(N)){
    N = length(dat$y)
    cat("N is not specified, defaulting to the number of observed values of y\n")
  }
  
  with(dat, 
       list(y=y, 
            X=as.matrix(select(dat, expnm)), 
            dx=dx, 
            p=p,
            N=N)
  )
}

convergence_check <- function(fit, ...){
#' @title lorem ipsum
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param fit ipsum
#' @param ... ipsum
#' @export
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
                         p = NULL,
                         dx = NULL,
                         N = NULL,
                         ...
                ){
#' @title lorem ipsum
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param i ipsum
#' @param raw ipsum
#' @param fl ipsum
#' @param s.code name of a character string with a stan model
#' @param s.file character name of a presepecified stan model ("wwbd_simtemplate_20190205", "wwbd_simlogistic_TransParameter_0211201920190211")
#' @param iter ipsum
#' @param warmup ipsum
#' @param p ipsum
#' @param dx ipsum
#' @param N ipsum
#' @param ... ipsum
#' @export
#' @import rstan
#' @examples
#' runif(1)
  #require(rstan)
  # do single analysis of data
  #stan model
  if(is.null(s.code) & is.null(s.file)) {
    cat("no stan model given, defaulting to horseshoe prior")
  s.code <- '
   // example with horseshoe prior
   data{
    int<lower=0> N;
    int<lower=0> p;
    int<lower=0> dx;
    vector[dx] X[N]; // functions like Nx5 matrix (but each column is real[])
    int y[N];
   }
   transformed data{
    //real meany;
    //vector[N] ycen;
    matrix[N,p] D;

    //meany = mean(y);
    //ycen = y-meany;
   // todo: transform y to center!
    for(c in 1:dx){
     D[,c] = to_vector(X[,c]);
     //D[,c] = X[,c];
    }
    for(c in (dx+1):p){
     D[,c] = to_vector(X[,c-5]) .* to_vector(X[,c-5]);
     //D[,c] = X[,c-5] .* X[,c-5];
    }
   }
   parameters{
    vector<lower=0>[p] lambda; // local shrinkage
    vector[p] beta;
    real<lower=0> sigma;
    real b0; // given uniform prior
    real<lower=0> tau; // global shrinkage
    //real<lower=0> sig; // global shrinkage hyperprior
   }
   transformed parameters{}
   model{
    {
     vector[N] mu;
     lambda ~ cauchy(0,1.0); // local shrinkage
     tau ~ cauchy(0,1.0); // global shrinkage
     target += -log(sigma*sigma); // jeffreys prior
     beta ~ normal(0,lambda * tau * sigma);
     mu = b0 + D * beta;
     y ~ bernoulli_logit(mu);
    }
   }
   generated quantities{
   real rd;
    {
     vector[N] r1;
     vector[N] r0;
     matrix[N,p] D1;
     matrix[N,p] D0;
      D1=D;
      D0=D;
     for(i in 1:N){
        // this is intervention to set exposure 1 to 1.0 versus 0.0 (i.e. both main effect
        // and the self interaction term go to 1.0
        D1[i,1] = 1.0;
        D1[i,6] = 1.0;
        D0[i,1] = 0.0;
        D0[i,6] = 0.0;
      }
    
      r1 =  inv_logit(b0 + D1 * beta);
      r0 =  inv_logit(b0 + D0 * beta);
      rd = mean(r1)-mean(r0);
    }
   }
'
  }
  sdat = data_reader(raw, i, p=p, dx=dx, N=N)
  # note: keep warming up even for a while after apparent convergence: helps improve efficiency
  #  of adaptive algorithm
  if(!is.null(s.file)){
     res = stan(file = system.file("stan", paste0(s.file,".stan", package = "wellwise")), 
             data = sdat, 
             chains = 4, 
             sample_file=fl, 
             iter=iter, 
             warmup=warmup, 
             ...)
  } else if(!is.null(s.code)){
         res = stan(model_code = s.code, 
             data = sdat, 
             chains = 4, 
             sample_file=fl, 
             iter=iter, 
             warmup=warmup, 
             ...)
  }
  #step to include here: automated monitoring for convergence and continued sampling if not converged
  #class(res) <- 'bgfsimmod'
  # inherets 'stanfit' class
  res
}

analysis_wrapper <- function(simiters, 
                             rawdata, 
                             dir="~/temp/", 
                             root, 
                             debug=FALSE, 
                             verbose=FALSE,
                             ...
                             ){
#' @title lorem ipsum
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param simiters ipsum
#' @param rawdata ipsum
#' @param dir ipsum
#' @param root ipsum
#' @param debug ipsum
#' @param verbose ipsum
#' @param ... ipsum
#' @export
#' @examples
#' runif(1)

  call = match.call()
  # perform multiple analyses
  if(length(simiters)==1) {
    sq = 1:simiters
  } else{
    sq = simiters
  }
  #mf = match(c("fl"), names(call))
  outfile = paste0(dir, root, "_res.csv")
  if(debug){
    outfile = "samples.csv"
    cat(paste0("Outputting samples from stan to ", outfile))
  }
  if(verbose) cat(paste0("Analyzing data ", length(sq), " times\n"))
  if(verbose) cat(paste0("R output can be seen at ", paste0(dir, root, "_rmsg.txt"), "\n"))
  res = list(1:length(sq))
  j=1
  sink(paste0(dir, root, "rmsg.txt"), split = FALSE, type = c("output", "message"))
  for(i in sq){
    cat(".")
    res[[j]] = data_analyst(i, rawdata, fl=outfile, ...)
    j=j+1
  }
  sink()
  cat("\n")
  class(res) <- 'bgfsimmodlist'
  res
}

summary.bgfsimmodlist <- function(
  object=NULL,
  ...
  ){
#' @name summary
#' @title lorem ipsum
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param object ipsum
#' @param ... ipsum
#' @import rstan
#' @method summary bgfsimmodlist
#' @importFrom stats sd quantile
#' @export
#' @return a 'bgfsimres' object
#' @examples
#' runif(1)
  # summarize
  # to do
  # if have values stored from analysis_wrapper, then use those
  # else read csv files
  numf = length(object)
  postmeans = list(numf)
  divergents = numeric(numf)
  for(r in 1:numf){
    tt = get_sampler_params(object[[r]])
    divergents[r] = sum(c(lapply(tt, function(x) sum(x[,"divergent__"]))>0))

    sum = summary(object[[r]])$summary
    idx = grep('rd', rownames(sum))
    postmeans[[r]] = summary(object[[r]])$summary[idx,1]
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
#' @title lorem ipsum
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
#' @title lorem ipsum
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
