#' Apply to HiC Experiment
#' 
#' @param X HiC experiment.
#' @param f Function to apply. 
#' @param ... Additional arguments for f.
#' @import doMC
#' @export

expApply = function(X,f,cores,...){
  # X : HiC experiment
  # f : Function mapping df to df applied to chromosome level data
  # Test : Testing ? 
  # ... : Additional arguments for f
  if(!missing(cores)){registerDoMC(cores=cores);}
  # Chromosomes
  chrs = names(X);
  n = length(chrs);
  
  # Loop over chromosomes
  Results = foreach(i=1:n) %dopar% {
    # Subset
    key = chrs[i];
    Sub = X[[key]];
    Out = f(Sub,...)
    return(Out);
  }
  names(Results) = names(X);
  # Output
  return(Results);
}

#' Expand HiC Experiment
#' 
#' Expands a collapsed HiC experiment into a list with 
#' one data frame per chromosome.
#' @param X HiC experiment.
#' @import foreach
#' @export

expandExp = function(X){
  # Chromosomes
  names = unique(X$Chr2);
  nchr = length(names);
  # Loop over chromosomes
  Out = foreach(i=1:nchr) %do% {
    Sub = X[X$Chr2==names[i],];
    return(Sub);
  }
  # Format
  names(Out) = names;
  return(Out);
}