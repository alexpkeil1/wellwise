
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
mkints <- function(terms, intvars){
#' @title lorem ipsum
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param terms ipsum
#' @param intvars ipsum
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




#create a model based on generic terms (mkmodidx probably better)
mkmod <- function(coef, terms, start=1){
#' @title lorem ipsum
#' 
#' @description lorem ipsum
#' 
#' @details lorem ipsum
#' @param coef ipsum
#' @param terms ipsum
#' @param start ipsum
#' @export
#' @examples
#' runif(1)
 mod <- ifelse(start==1, paste0(coef,"0"), "")
 for(i in (start-1)+1:length(terms)){
   mod <- paste0(mod, " + ", coef, "[", i, "]*", terms[i-(start-1)])
 }
cat( "\n", mod, "\n\n")
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
   mod <- gsub(paste0(var, "([\\[ ])"), paste0(sub, "\\1"), paste0(mod, " "))
 }
 trim(mod)
}



