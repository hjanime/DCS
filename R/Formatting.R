#' @include ClassDefinitions.R
NULL

# Purpose : Functions to convert sparse HiC objects into DCS objects

#' Sparse Matrix to DCS Chromosome
#' 
#' @param A Sparse matrix of class dtCMatrix
#' @param chr Chromosome
#' @param smp Sample name
#' @param group Group name
#' @return A contact matrix of class \code{DCSchr}.
#' 
#' @importClassesFrom Matrix dtCMatrix
#' @importFrom Matrix summary t
#' @importFrom methods new

sparseMat2DCS = function(A,chr,smp,group){
  # Check class
  if(!(class(A) %in% c("dtCMatrix","dtTMatrix"))){stop("Expecting a sparse matrix.")};  
  # Convert to data frame
  B = summary(t(A));
  # Coord frame
  Coord = B[,c("i","j")];
  Coord$d = (Coord$j-Coord$i);
  Coord = as.matrix(Coord);
  # Count frame
  Counts = data.matrix(B[,3,drop=F]);
  colnames(Counts) = smp;
  # Output structure
  Out = new(Class="DCSchr",chr=as.character(chr),Coord=Coord,Counts=Counts,
            groups=group,samples=as.character(smp),Stats=data.frame());
  return(Out);
}

#' Convert List of Sparse Matrices to a DCS Experiment
#' 
#' @param L List of sparse matrices
#' @param chr Vector of chromosome names, one per sparse matrix
#' @param smp Optional sample name
#' @param group Optional group assignment
#' @return An object of class \code{DCSexp}.
#' 
#' @importFrom foreach "%do%" foreach
#' @importFrom methods new
#' @export 
#' 
#' @examples 
#' # Auxin 
#' A1 = DCS::Auxin_rep1;
#' A2 = DCS::Auxin_rep2;
#' # Combine sparse matrices for individual chromosomes within a sample
#' Auxin1 = DCS::createDCSsmp(L=A1,chr=c("1","2"),smp="Aux.R1",group=0);
#' Auxin2 = DCS::createDCSsmp(L=A2,chr=c("1","2"),smp="Aux.R2",group=0);

createDCSsmp = function(L,chr,smp,group){
  nc = length(chr);
  if(nc!=length(L)){stop("Must supply as many chromosome names as sparse matrices.")};
  if(length(unique(chr))!=length(L)){stop("Chromosome names must be unique.")};
  if(!missing(smp)){smp = as.character(smp)} else {smp = "Sample"};
  if(!missing(group)){groups = as.factor(group);} else {groups = factor(0)};
  # Create list of DCS chromosomes
  i = NULL;
  Lout = foreach(i=1:nc) %do% {
    return(sparseMat2DCS(A=L[[i]],chr=chr[i],smp=smp,group=groups));
  }
  names(Lout) = chr;
  Out = new("DCSexp",chrs=chr,Data=Lout,groups=groups,samples=smp);
  return(Out);
}

#' Combine Multiple Samples into a DCS Experiment
#' 
#' @param L List of DCSexp objects corresponding to different samples
#' @param smps Sample names, required if not already assigned to component samples
#' @param groups Optional group assignments
#' @return An object of class \code{DCSexp}.
#' 
#' @importFrom foreach "%do%" foreach
#' @importFrom methods new
#' @export
#' 
#' @examples
#' # Auxin 
#' A1 = DCS::Auxin_rep1;
#' A2 = DCS::Auxin_rep2;
#' # Control
#' C1 = DCS::Control_rep1;
#' C2 = DCS::Control_rep2;
#' # Combine sparse matrices for individual chromosomes within a sample
#' Auxin1 = DCS::createDCSsmp(L=A1,chr=c("1","2"),smp="Aux.R1",group=0);
#' Auxin2 = DCS::createDCSsmp(L=A2,chr=c("1","2"),smp="Aux.R2",group=0);
#' Ctrl1 = DCS::createDCSsmp(L=C1,chr=c("1","2"),smp="Ctrl.R1",group=1);
#' Ctrl2 = DCS::createDCSsmp(L=C2,chr=c("1","2"),smp="Ctrl.R2",group=1);
#' # Combine samples into an experiment
#' Exp = DCS::createDCSexp(L=list(Auxin1,Auxin2,Ctrl1,Ctrl2));

createDCSexp = function(L,smps,groups){
  # Number of samples
  N = length(L);
  # Obtain sample names if not provided
  if(missing(smps)){
    aux = function(x){x@samples};
    smps = sort(unlist(lapply(L,aux)));
  }
  # Check for as many unique sample names as samples
  if(length(unique(smps))!=N){stop("Must supply as many unique sample names as samples.")}
  # Obtain sample groups if not provided
  if(missing(groups)){
    aux = function(x){x@groups};
    groups = sort(unlist(lapply(L,aux)));
  }
  # Ensure all samples have the same chromosomes
  aux = function(x){x@chrs};
  chrs = lapply(L,aux);
  chrs.equal = all.equal(Reduce(chrs,f=intersect),chrs[[1]]);
  if(!(chrs.equal==T)){stop("All samples must contain the sample chromosomes.")};
  # Number of chromosomes
  chrs = sort(chrs[[1]]);
  nc = length(chrs);
  # Loop over chromosomes
  i = j = NULL;
  Lout = foreach(i=1:nc) %do% {
      # Loop over samples
      Sub = foreach(j=1:N) %do% {
        A = L[[j]]@Data[[chrs[i]]];
      }
      # Merge samples
      Out = Reduce(Sub,f=merge.DCSchr);
  }
  # Output
  names(Lout) = chrs;
  Out = new("DCSexp",chrs=chrs,Data=Lout,groups=groups,samples=smps);
  return(Out);
}