
#helper function
trim <- function (x){
#' @title lorem ipsum
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param x ipsum
#' @export
#' @examples
#' runif(1)
 gsub("^\\s+|\\s+$", "", x)
}

#helper function
sortstring <- function(s){
#' @title lorem ipsum
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param s ipsum
#' @export
#' @examples
#' runif(1)
 paste(sort(unlist(strsplit(s, "\\.\\*"))), collapse = ".*")
}  
  
#create interaction terms
mkints <- function(terms, intvars, binvars){
#' @title lorem ipsum
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param terms ipsum
#' @param intvars ipsum
#' @param binvars binary variables (character vector) which are excluded to avoid meaningless self-interactions
#' @export
#' @examples
#' runif(1)
 out <- ""
 for(i in 1:length(intvars)){
  newout <- paste(terms, intvars[i], sep=".*")
  for(j in 1:length(newout)){
    newout[j] <- sortstring(newout[j])
  }
  out <- unique(c(out, newout))
 }
 if(!is.null(binvars)){
  excl = paste(binvars, binvars, sep='.*')
  out = setdiff(out, excl)
 }
 trim(paste(out, collapse=" "))
}

#trim some self-interaction variables (say x*x if x is binary)
trimints <- function(trimvars, terms){
#' @title lorem ipsum
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param trimvars ipsum
#' @param terms ipsum
#' @export
#' @examples
#' runif(1)
  newterms=terms
  for(i in 1:length(trimvars)){
   trimterm <- paste0(trimvars[i], "\\.\\*",trimvars[i])
   newterms <- gsub(trimterm, "", newterms)
 }
 newterms[newterms!=""]
}


make_fullrmod <- function(terms, intvars, binvars=NULL){
#' @title make a model to use in model.matrix
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param terms ipsum
#' @param intvars ipsum
#' @param binvars ipsum
#' @export
#' @importFrom stats as.formula
#' @examples
#' runif(1)
 mod <- paste(terms, collapse = " + ")
 mod2 <- mkints(terms=terms, intvars=terms, binvars=NULL)
 f = paste("~ -1 + ", mod, " + ", gsub(' ', ' + ', gsub('\\.', '', mod2)), collapse=' + ')
 f = gsub("([A-Za-z]+)\\*([A-Za-z]+)", "I(\\1*\\2)", f)
 as.formula(f)
}



#create a model based on generic terms (mkmodidx probably better)
mkmod <- function(coef, terms, start=1, catmod=FALSE){
#' @title lorem ipsum
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param coef ipsum
#' @param terms ipsum
#' @param start ipsum
#' @param catmod ipsum
#' @export
#' @examples
#' runif(1)
 mod <- ifelse(start==1, paste0(coef,"0"), "")
 for(i in (start-1)+1:length(terms)){
   mod <- paste0(mod, " + ", coef, "[", i, "]*", terms[i-(start-1)])
 }
if(catmod) cat( "\n", mod, "\n\n")
return(mod)
}

mkdesignR <- function(terms, catmod=FALSE){
#' @title lorem ipsum
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param terms ipsum
#' @param catmod ipsum
#' @export
#' @examples
#' runif(1)
 mod <- paste0("rep(1,n)")
 for(i in 1:length(terms)){
   mod <- paste0(mod, ", ", gsub("\\.","",terms[i]))
 }
if(catmod) cat( "\n", mod, "\n\n")
return(mod)
}

#create a model based on generic terms, including interactions, and output model as a string
mkmodidx <- function(coef, terms, start=1, index="i", print=FALSE, indexed=FALSE){
#' @title lorem ipsum
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param coef ipsum
#' @param terms ipsum
#' @param start ipsum
#' @param index ipsum
#' @param print ipsum
#' @param indexed ipsum
#' @export
#' @examples
#' runif(1)
  mod <- ifelse(start==1, paste0(coef,"0"), "")
  for(i in (start-1)+1:length(terms)){
    mod <- paste0(mod, " + ", coef, "[", i, "]*", terms[i-(start-1)])
  }
  mod2 <- gsub("([_a-z]+[a-z0-9_]*)", paste0("\\1\\[", index, "\\]"), mod)
  mod2 <- gsub(paste0("(\\[", index, "\\])(\\[[0-9]+\\])"), "\\2", mod2)
  mod2 <- gsub(paste0("(",coef,")(\\[", index, "\\])(0)"), "\\1\\3", mod2)
  if(print){
   cat("\n", "Vectorized")
   cat( "\n", mod, "\n\n")
   cat("\n", "Indexed")
   cat( "\n", mod2, "\n\n")
  }
  if(indexed) return(gsub("\\.(\\*)", "\\1",mod2)) 
  else return(gsub("(\\.)(\\*)", " (\\1) (\\2) ",mod))
}

mkintervention <- function(mod, vars, subs){
#' @title lorem ipsum
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param mod ipsum
#' @param vars ipsum
#' @param subs ipsum
#' @export
#' @examples
#' runif(1)
 for(i in 1:length(vars)){
  var <- vars[i]
  sub <- subs[i]
   mod <- gsub(paste0(var, "([\\[ -])"), paste0(sub, "\\1"), paste0(mod, " "))
 }
 trim(mod)
}


mkinterventionR <- function(mod, vars, subs){
#' @title lorem ipsum
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param mod ipsum
#' @param vars ipsum
#' @param subs ipsum
#' @export
#' @examples
#' runif(1)
 for(i in 1:length(vars)){
  var <- vars[i]
  sub <- subs[i]
   mod <- gsub(paste0(var), paste0(sub), paste0(mod, ""))
 }
 trim(mod)
}


stan_basic <- function(x=c('x', 'z'),
                       z=c('bmi'),
                       y='y',
                       binvars = NULL,
                       matx="X",
                       standardizex=TRUE,
                       binary=TRUE,
                       vectorized=FALSE,
                       xintv=rbind(c(.99, 0),c(0, .99),c(.99, .99))){
#' @title make a basic stan model
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param x intervenable exposures (character vector)
#' @param z covariates (character vector)
#' @param y outcome
#' @param binvars = non-outcome variables that are binary (character vector)
#' @param matx optional, name of matrix with intervenable exposures
#' @param standardizex logical, should x be standardized?
#' @param binary logical, is outcome binary?
#' @param vectorized logical, should model be expressed in vector notation when possible?
#' @param xintv matrix with ncol = number of intervenable exposures, nrow = number of interventions. Each value is on [0,1] and represents the proportional decrease in the value of x upon hypothetical intervention
#' @importFrom stringr str_split str_wrap
#' @export
#' @examples
#' # library(rstan)
#'  dgm <- function(N=100, trueRD=0.2){
#'    x1 = rbinom(N, 1, 0.5)
#'    py00 = runif(N)*0.1 + 0.4
#'    l2 = rbinom(N, 1, 1/(1+exp(-1 + x1 + py00)))
#'    x2 = rbinom(N, 1, 1/(1+exp(-1 + x1 + l2)))
#'    py = py00 + trueRD*((x1 + x2)/2) #true risk difference per unit exposure;
#'    y = rbinom(N, 1, py)
#'    data.frame(x1, l2, x2, y)
#'  }
#'  dat = as.list(dgm(100))
#'  dat$N = 100
#'  dat$p = 5
#'  
#'  source("~/Epiprojects/wellwater/sims/code/make_stan_terms.R")
#'  
#'  mod = stan_basic(x=c('x1', 'x2'), z = 'l2', y='y', 
#'    binvars=c('x1', 'x2', 'l2'), xintv = rbind(c(1,0), c(0,1),c(1,1)), 
#'    vectorized = TRUE, binary=TRUE, matx = NULL)
#'  cat(mod)
#' # usage in stan (or edit by hand)
#' # not run
#' # stan(model_code = mod, data = dat, chains=1, iter=100)
  # initialize each block  
  data <- 'data {\n  int<lower=0> N;\n  int<lower=0> p;'
  tdata <- 'transformed data {'
  parms <- 'parameters {'
  tparms <- 'transformed parameters {'
  model <-  'model {'
  gquant <- 'generated quantities {'
  
  # data
  if(binary){
    data = paste0(data, "\n  int<lower=0, upper=1> ", y, "[N];")
  }
  if(!binary){
    data = paste0(data, "\n  real ", y, "[N];")
  }
  if(!vectorized){
    #if(is.null(matx)) for(n in x) data = (paste0(data, "\n  real ", n, "[N];"))
    #for(n in z) data = (paste0(data, "\n  real ", z, "[N];"))
    if(is.null(matx)) for(n in x) data = (paste0(data, "\n  vector[N] ", n, ";"))
    for(n in z) data = (paste0(data, "\n  vector[N] ", n, ";"))
  }
  if(vectorized){
    if(is.null(matx)) for(n in x) data = (paste0(data, "\n  vector[N] ", n, ";"))
    for(n in z) data = (paste0(data, "\n  vector[N] ", n, ";"))
  }
  if(!is.null(matx)){
    data = paste0(data, "\n  matrix[N,", length(x),"] ", matx, ";")
  }
  # transformed data
  if(!is.null(matx)){
    j = 1
    for(n in x) {
      tdata = paste0(tdata, "\n  vector[N] ", n, " = col(", matx, ",",j,");")
      j=j+1
    }
  }
  
   # standardized non-binary exposures
   if(standardizex){
     for(ix in 1:length(x)){
      if(!(x[ix] %in% binvars)){
       tdata = paste0(tdata, paste0("\n  vector[N] cen_", x[ix], ";"))
       tdata = paste0(tdata, paste0("\n  real m", x[ix], " = mean(", x[ix], ");"))
       tdata = paste0(tdata, paste0("\n  real s", x[ix], " = sd(", x[ix], ");"))
      }
     }
     for(ix in 1:length(x)){
      if(!(x[ix] %in% binvars)){
       tdata = paste0(tdata, paste0("\n  cen_", x[ix], " = (", x[ix] ,"-m",x[ix] ,") ./ s",x[ix],";"))
      }
     }     
     ox = x # original x
     xs = x
     for(ix in 1:length(x)){
      if(!(x[ix] %in% binvars)){
        xs[ix] = paste0("((", x[ix] ,"-m",x[ix] ,") ./ s",x[ix],")") # x with long hand standardization
        x[ix] = paste0('cen_', x[ix])  # standardized x
      }
     }
   }
   if(!standardizex){
     ox=x
     xs=x
   }   
  # parameters
    parms = paste0(parms, '\n  vector[p] b_;')
    parms = paste0(parms, '\n  real beta0;') 
    if(!binary) parms = paste0(parms, '\n  real<lower=0> sigma;') 
  # transformed parameters
    tparms = paste0(tparms, '\n  vector[p] beta;')
    tparms = paste0(tparms, '\n  beta = b_ * 10.0;')
  
  # model
  model = paste0(model, '\n  {\n  vector[N] mu;;')
  intx = str_split(mkints(terms=c(x,z), intvars=c(x,z), binvars=binvars), " ")[[1]]
  model = paste0(model, '\n  beta0 ~ normal(0, 10);')
  model = paste0(model, '\n  b_ ~ normal(0, 1);')
  if(!binary) model = paste0(model, '\n  sigma ~ cauchy(0, 1);//// half cauchy')
  if(!vectorized){
    model = paste0(model, '\n   for(i in 1:N){')
    mucode <- mkmodidx(coef='beta', terms=c(c(x,z), intx), start=1, print = FALSE, indexed = TRUE)
    mucode <- gsub('beta0[i]', 'beta0', mucode, fixed=TRUE)
    model = paste0(model, '\n    mu[i] = ', str_wrap(mucode, 80, exdent=12), ';')
    if(binary) model = paste0(model, '\n   y ~ bernoulli_logit(mu[i]);')
    if(!binary) model = paste0(model, '\n   y ~ normal(mu[i], sigma);')
    model = paste0(model, '\n   }//i')
  }
  if(vectorized){
    m1 = mkmod(coef='beta', terms=c(c(x,z), intx), start=1) # better
    mucode = gsub("\\.\\*", " .*", m1)
    model = paste0(model, '\n  mu = ', str_wrap(mucode, 80, exdent=11), ';')
    if(binary) model = paste0(model, '\n   y ~ bernoulli_logit(mu);')
    if(!binary) model = paste0(model, '\n   y ~ normal(mu, sigma);')
  }
  model = paste0(model, '\n    }//local')


  # generated quantities
  if(is.null(xintv)) xintv = rbind(rep(0, length(x)))
  gquant = paste0(gquant, '\n  vector[', nrow(xintv),'] rd;\n {')
  gquant = paste0(gquant, '\n  vector[N] r1;\n  vector[N] r0;')
  
  gquant = paste0(gquant, '\n  matrix[', nrow(xintv), ',', ncol(xintv),'] intprop = [\n')
  gquant = paste0(gquant, paste('      [', apply(xintv, 1, function(x) paste(x, collapse=',')), ']', collapse=',\n'))
  gquant = paste0(gquant, '];\n')
  # use standardized versions
  if(standardizex){
    for(j in 1:length(x)) mucode = gsub(x[j], xs[j], mucode, fixed=TRUE) # replace standardized with long hand standardized
    if(!vectorized){
      for(j in 1:length(x)) mucode = gsub(paste0(xs[j], '[i]'), xs[j], mucode, fixed=TRUE) # replace standardized with long hand standardized
      for(j in 1:length(x)) mucode = gsub(paste0(ox[j], '-'), paste0(ox[j], '[i]-'), mucode, fixed=TRUE) # replace standardized with long hand standardized
      #binvars
      for(xr in intersect(binvars, x))
        mucode = gsub(xr, paste0(xr, '[i]'), mucode, fixed=TRUE) # replace standardized with long hand standardized
    }
  }
  if(vectorized){
    nc = mkintervention(mucode, vars=ox, subs=ox) 
    gquant = paste0(gquant, '\n  r0 = ', str_wrap(nc, 80, exdent=11), ';')
    gquant = paste0(gquant, '\n   for(j in 1:',nrow(xintv),'){')
    subs = sapply(1:ncol(xintv), function(i) paste0('((1-intprop[j,', i, '])*',ox[i],')'))
    ints = mkintervention(mucode, vars=ox, subs=c(subs)) 
    gquant = paste0(gquant, '\n    r1 = ', str_wrap(ints, 80, exdent=11), ';')
    if(binary) gquant = paste0(gquant, '\n    rd[j] = mean(inv_logit(r1))-mean(inv_logit(r0));')
    if(!binary) gquant = paste0(gquant, '\n    rd[j] = mean(r1)-mean(r0);')
    gquant = paste0(gquant, '\n  }//j')
  }
  if(!vectorized){
    gquant = paste0(gquant, '\n  for(j in 1:',nrow(xintv),'){')
    gquant = paste0(gquant, '\n    for(i in 1:N){')
    nc = mkintervention(mucode, vars=ox, subs=ox) 
    gquant = paste0(gquant, '\n      r0[i] = ', str_wrap(nc, 80, exdent=12), ';')

    subs = sapply(1:ncol(xintv), function(i) paste0('(1-intprop[j,', i, '])*',ox[i],''))
    ints = mkintervention(mucode, vars=ox, subs=c(subs)) 
    gquant = paste0(gquant, '\n      r1[i] = ', str_wrap(ints, 80, exdent=12), ';')
    gquant = paste0(gquant, '\n    }//i')
    if(binary) gquant = paste0(gquant, '\n    rd[j] = mean(inv_logit(r1))-mean(inv_logit(r0));')
    if(!binary) gquant = paste0(gquant, '\n    rd[j] = mean(r1)-mean(r0);')
    gquant = paste0(gquant, '\n  }//j')
  } 

  # end with a bracket  
  data <- paste0(data, '\n}// end data')
  tdata <- paste0(tdata, '\n}// end transformed data')
  parms <- paste0(parms, '\n}// end parameters')
  tparms <- paste0(tparms, '\n}// end transformed parameters')
  model <- paste0(model, '\n}// end model')
  gquant <- paste0(gquant, '\n }//local\n}// end stan code\n\n')
  paste(data, tdata, parms, tparms, model, gquant, sep = "\n")
}



jags_basic <- function(x=c('x', 'z'),
                       z=c('bmi'),
                       y='y',
                       binvars = NULL,
                       matx="X",
                       standardizex=TRUE,
                       binary=TRUE,
                       #vectorized=FALSE,
                       xintv=rbind(c(.99, 0),c(0, .99),c(.99, .99))){
#' @title make a basic jags model
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param x intervenable exposures (character vector)
#' @param z covariates (character vector)
#' @param y outcome
#' @param binvars = non-outcome variables that are binary (character vector)
#' @param matx optional, name of matrix with intervenable exposures
#' @param standardizex logical, should x be standardized?
#' @param binary logical, is outcome binary?
#' @param xintv matrix with ncol = number of intervenable exposures, nrow = number of interventions. Each value is on [0,1] and represents the proportional decrease in the value of x upon hypothetical intervention
#' @importFrom stringr str_split str_wrap
#' @export
#' @examples
#' # library(rjags)
#'  dgm <- function(N=100, trueRD=0.2){
#'    x1 = rbinom(N, 1, 0.5)
#'    py00 = runif(N)*0.1 + 0.4
#'    l2 = rbinom(N, 1, 1/(1+exp(-1 + x1 + py00)))
#'    x2 = rbinom(N, 1, 1/(1+exp(-1 + x1 + l2)))
#'    py = py00 + trueRD*((x1 + x2)/2) #true risk difference per unit exposure;
#'    y = rbinom(N, 1, py)
#'    data.frame(x1, l2, x2, y)
#'  }
#'  dat = as.list(dgm(100))
#'  dat$N = 100
#'  dat$p = 5
#'  
#'  source("~/Epiprojects/wellwater/sims/code/make_stan_terms.R")
#'  
#'  mod = jags_basic(x=c('x1', 'x2'), z = 'l2', y='y', 
#'    binvars=c('x1', 'x2', 'l2'), xintv = rbind(c(1,0), 
#'    c(0,1),c(1,1)), binary=TRUE, matx = NULL)
#'  cat(mod)
#' # usage in jags (or edit by hand)
#' # not run
#' # tf = tempfile()
#' # cat(mod, file=tf)
#' # jags.model(file = tf, data = dat, n.chains=1)
  data <- 'data {\n  dummy <- 1'
  model <-  'model {'
  
  if(!is.null(matx)){
    j = 1
    for(n in x) {
      data = paste0(data, "\n  ", n, " <- ", matx, "[,",j,"]")
      j=j+1
    }
  }

  
  intx = str_split(mkints(terms=c(x,z), intvars=c(x,z), binvars=binvars), " ")[[1]]
  #likelihood
  if(!binary) model = paste0(model, '\n  sigma ~ cauchy(0, 1);# half cauchy')
  model = paste0(model, '\n  for(i in 1:N){')
  mucode <- mkmodidx(coef='beta', terms=c(c(x,z), intx), start=1, print = FALSE, indexed = TRUE)
  mucode <- gsub('beta0[i]', 'beta0', mucode, fixed=TRUE)
  if(binary) {
    model = paste0(model, '\n    y[i] ~ dbern(mu[i]);')
    model = paste0(model, '\n    mu[i] <- ilogit(', str_wrap(mucode, 80, exdent=12), ');')
  }
  if(!binary){
    model = paste0(model, '\n    y[i] ~ dnorm(mu[i], sigma);')
    model = paste0(model, '\n    mu[i] <- ', str_wrap(mucode, 80, exdent=12), ';')
  }
  # interventions
   if(!standardizex){
     ox=x
     xs=x
   }   
   if(standardizex){
     # do nothing, for now
     ox=x
     xs=x
   }   

  if(is.null(xintv)) xintv = rbind(rep(0, length(x)))
  for(ridx in 1:nrow(xintv)){
    #for(cidx in 1:ncol(xintv)){
    #  data = paste0(data,
    #                '\n  intprop[', ridx, ',',cidx, '] = ',xintv[ridx,cidx] ,'')
    #}
    data = paste0(data,'\n  intprop[', ridx, ',1:',ncol(xintv),'] = c(',paste0(xintv[ridx,], collapse=',') ,')')    
  }
  subs = sapply(1:ncol(xintv), function(i) paste0('(1-intprop[j,', i, '])*',ox[i],''))
  ints = mkintervention(mucode, vars=ox, subs=c(subs)) 
  #intervention code
   model = paste0(model, '\n    for(j in 1:',nrow(xintv), '){' )
   if(binary)  model = paste0(model, '\n      r1[i,j] = ilogit(', str_wrap(ints, 80, exdent=12), ');')
   if(!binary) model = paste0(model, '\n      r1[i,j] = ', str_wrap(ints, 80, exdent=12), ';')
   model = paste0(model, '\n    }#j')
  
  model = paste0(model, '\n  }#i')
  model = paste0(model, '\n  for(j in 1:',nrow(xintv), '){' )
  model = paste0(model, '\n    rd[j]= mean(r1[,j])-mean(mu)')
  model = paste0(model, '\n  }#j')
  #priors
  model = paste0(model, '\n  beta0 ~ dnorm(0, 10);')
  model = paste0(model, '\n  for(m in 1:p){')
  model = paste0(model, '\n    beta[m] ~ dnorm(0, 1);')
  model = paste0(model, '\n  }#m')
  
  data <- paste0(data, '\n}# end data')
  model <- paste0(model, '\n# end model\n}')
  paste(data, model, sep = "\n")
}



gibbs_basic <- function(x=c('x', 'z'),
                       z=c('bmi'),
                       y='y',
                       binvars = NULL,
                       matx="X",
                       standardizex=TRUE,
                       binary=TRUE,
                       #vectorized=FALSE,
                       xintv=rbind(c(.99, 0),c(0, .99),c(.99, .99))){
#' @title make a basic jags model
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param x intervenable exposures (character vector)
#' @param z covariates (character vector)
#' @param y outcome
#' @param binvars = non-outcome variables that are binary (character vector)
#' @param matx optional, name of matrix with intervenable exposures
#' @param standardizex logical, should x be standardized?
#' @param binary logical, is outcome binary?
#' @param xintv matrix with ncol = number of intervenable exposures, nrow = number of interventions. Each value is on [0,1] and represents the proportional decrease in the value of x upon hypothetical intervention
#' @importFrom stringr str_split str_wrap
#' @export
#' @examples
#' # library(rjags)
#'  dgm <- function(N=100, trueRD=0.2){
#'    x1 = rbinom(N, 1, 0.5)
#'    py00 = runif(N)*0.1 + 0.4
#'    l2 = rbinom(N, 1, 1/(1+exp(-1 + x1 + py00)))
#'    x2 = rbinom(N, 1, 1/(1+exp(-1 + x1 + l2)))
#'    py = py00 + trueRD*((x1 + x2)/2) #true risk difference per unit exposure;
#'    y = rbinom(N, 1, py)
#'    data.frame(x1, l2, x2, y)
#'  }
#'  dat = as.list(dgm(100))
#'  dat$N = 100
#'  dat$p = 5
#'  
#'  source("~/Epiprojects/wellwater/sims/code/make_stan_terms.R")
#'  
#'  mod = jags_basic(x=c('x1', 'x2'), z = 'l2', y='y', 
#'    binvars=c('x1', 'x2', 'l2'), xintv = rbind(c(1,0), c(0,1),c(1,1)), 
#'    binary=TRUE, matx = NULL)
#'  cat(mod)
#' # usage in jags (or edit by hand)
#' # not run
#' # tf = tempfile()
#' # cat(mod, file=tf)
#' # jags.model(file = tf, data = dat, n.chains=1)
  data <- ''
  model <-  ''
  
  if(!is.null(matx)){
    j = 1
    for(n in x) {
      data = paste0(data, "\n  ", n, " <- ", matx, "[,",j,"]")
      j=j+1
    }
  }

  intx = str_split(mkints(terms=c(x,z), intvars=c(x,z), binvars=binvars), " ")[[1]]
  intx = gsub("//.", "", intx)
  #likelihood
  mucode <- mkdesignR(terms=c(c(x,z), intx))

  model = paste0(model, '\n    Xi = cbind(', str_wrap(mucode, 80, exdent=8), ');')
  # interventions
   if(!standardizex){
     ox=x
     xs=x
   }   
   if(standardizex){
     # do nothing, for now
     ox=x
     xs=x
   }   

  if(is.null(xintv)) xintv = rbind(rep(0, length(x)))
  for(ridx in 1:nrow(xintv)){
    data = paste0(data,'\n  intprop[', ridx, ',1:',ncol(xintv),'] = c(',paste0(xintv[ridx,], collapse=',') ,')')    
  }
  subs = sapply(1:ncol(xintv), function(i) paste0('(1-intprop[j,', i, '])*',ox[i],''))
  ints = mkinterventionR(mucode, vars=ox, subs=c(subs)) 
  #intervention code
   model = paste0(model, '\n    for(j in 1:',nrow(xintv), '){' )
   model = paste0(model, '\n      Xil[j] = cbind(', str_wrap(ints, 80, exdent=10), ');')
   model = paste0(model, '\n    }#j')
  
  model = paste0(model, '\n  for(j in 1:',nrow(xintv), '){' )
  model = paste0(model, '\n    rd[j]= mean(r1[,j])-mean(mu)')
  model = paste0(model, '\n  }#j')
  #priors

  data <- paste0(data, '\n# end data')
  model <- paste0(model, '\n# end model\n}')
  paste(data, model, sep = "\n")
}


############ usage ############