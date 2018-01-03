# Purpose: Define S4 classes used in DCS package

########################
# 1. Class definitions
########################

# Function to check DCSchr class validity
checkDCSchr = function(object){
  # Store errors
  errors = character();
  # 1. Structure of Coord matrix
  if(sum(colnames(object@Coord) != c("i","j","d"))>0){
    msg = "Coord matrix should have columns (i,j,d).";
    errors = c(errors,msg);
  }
  if(sum(object@Coord[,"i"]>object@Coord[,"j"])>0){
    msg = "Coord matrix requires i <= j."
    errors = c(errors,msg);
  }
  d = (object@Coord[,"j"]-object@Coord[,"i"]);
  if(sum(d != object@Coord[,"d"])>0){
    msg = "Distance d in Coord matrix does not equal (j-i)."
    errors = c(errors,msg);
  }
  # 2. Structure of Counts matrix
  if(sum(colnames(object@Counts) %in% object@samples) != length(object@samples)){
    msg = "Sample names do not correspond to columns of Count matrix."
    errors = c(errors,msg);
  }
  # Output
  if(length(errors)==0){return(TRUE)} else {return(errors)};
}

# Function to check DCSexp class validity
checkDCSexp = function(object){
  # Store errors
  errors = character();
  # 1. Check number of chromosomes corresponds to length of data list
  nc = length(object@chrs);
  nd = length(object@Data);
  if(nc!=nd){
    msg = "Mismatch between number of chromosomes and length of the data list."
    errors = c(errors,msg);
  }
  # 2. Check chromosome names correspond to elements of the data list
  if(sum(object@chrs %in% names(object@Data))!=nd){
    msg = "Chromosome names do not correspond with elements of the data list."
    errors = c(errors,msg);
  }
  # 3. Check type of elements in the data list
  if(sum(unlist(lapply(object@Data,FUN=class))=="DCSchr")!=nd){
    msg = "Elements of data list must have class DCSchr.";
    errors = c(errors,msg);
  }
  # Output
  if(length(errors)==0){return(TRUE)} else {return(errors)};
}

# Function to check FSchr class validity
checkFSchr = function(object){
  # Store errors
  errors = character();
  # 1. Check number of chromosomes corresponds to length of data list
  nc = length(object@groups);
  ns = length(object@Stats);
  if(nc!=ns){
    msg = "Mismatch between number of groups and length of the stats list."
    errors = c(errors,msg);
  }
  # 2. Check chromosome names correspond to elements of the data list
  if(sum(object@groups %in% names(object@Stats))!=ns){
    msg = "Group names do not correspond with elements of the stats list."
    errors = c(errors,msg);
  }
  # 3. Check type of elements in the data list
  # Output
  if(length(errors)==0){return(TRUE)} else {return(errors)};
}

# Function to check FSexp class validity
checkFSexp = function(object){
  # Store errors
  errors = character();
  # 1. Check number of chromosomes corresponds to length of data list
  nc = length(object@chrs);
  nd = length(object@Data);
  if(nc!=nd){
    msg = "Mismatch between number of chromosomes and length of the data list."
    errors = c(errors,msg);
  }
  # 2. Check chromosome names correspond to elements of the data list
  if(sum(object@chrs %in% names(object@Data))!=nd){
    msg = "Chromosome names do not correspond with elements of the data list."
    errors = c(errors,msg);
  }
  # 3. Check type of elements in the data list
  if(sum(unlist(lapply(object@Data,FUN=class))=="FSchr")!=nd){
    msg = "Elements of data list must have class FSchr.";
    errors = c(errors,msg);
  }
  # Output
  if(length(errors)==0){return(TRUE)} else {return(errors)};
}

# Class definitions

#' DCS Chromosome
#' 
#' DCS object for storing chromosome level data.
#' @slot chr Chromosome, character.
#' @slot Coord Matrix of pixel coordinates.
#' @slot Counts Matrix of pixel intensities.
#' @slot groups Group assignments, factor.
#' @slot samples Sample names, character.
#' @slot Stats Statistics available for each pixel, data.frame. 
#' @name DCSchr-class
#' @rdname DCSchr-class
#' @exportClass DCSchr
setClass(Class="DCSchr",representation=representation(chr="character",Coord="matrix",Counts="matrix",
                                 groups ="factor",samples="character",Stats="data.frame"),validity=checkDCSchr);

#' DCS Experiment
#' 
#' DCS object for storing experiment level data. 
#' @slot chrs Chromosomes in experiment, character.
#' @slot Data List of DCSchrs.
#' @slot groups Group assignments, factor.
#' @slot samples Sample names, character.
#' @name DCSexp-class
#' @rdname DCSexp-class
#' @exportClass DCSexp
setClass(Class="DCSexp",representation=representation(chrs="character",Data="list",groups="factor",samples="character"),
         validity=checkDCSexp);

#' FS Chromosome
#' 
#' Object for storing focal statistics (FS) at the chromosome level.
#' @slot chr Chromosome, character.
#' @slot foci Foci, or loci along the diagonal of the contact matrix, for which statistics are available.
#' @slot groups Group assignments, factor.
#' @slot Stats List of statistics, one per group. 
#' @name FSchr-class
#' @rdname FSchr-class
#' @exportClass FSchr
setClass(Class="FSchr",representation=representation(chr="character",foci="numeric",groups="factor",Stats="list"),validity=checkFSchr);

#' FS Experiment
#' 
#' Object for storing focal statistics (FS) at the experiment level.
#' @slot chrs Chromosomes, character.
#' @slot groups Group assignments, factor.
#' @slot Data List of FSchrs. 
#' @name FSexp-class
#' @rdname FSexp-class
#' @exportClass FSexp
setClass(Class="FSexp",representation=representation(chrs="character",groups="factor",Data="list"),validity=checkFSexp);

########################
# Apply Methods
########################

#' Apply to DCS Experiment
#' 
#' @param X A \code{DCSexp} object.
#' @param f Function to apply to each chromosome in \code{DCSexp@Data}. Should 
#'   accept objects of class \code{DCSchr}.
#' @param exp.out Is the output an object of class \code{DCSexp}? Default is 
#'   TRUE.
#' @param parallel Run in parallel? Must register parallel backend first. 
#'   Default is FALSE.
#' @param ... Additional arguments for f.
#' @return List resulting from application of function f to the chromosomes in
#'   \code{DCSexp@Data}. If \code{exp.out=T}, the list is formatted as a
#'   \code{DCSchr}.
#'   
#' @importFrom foreach "%dopar%" foreach registerDoSEQ
#' @export

apply.DCSexp = function(X,f,exp.out=T,parallel=F,...){
  # Chromosomes
  chrs = X@chrs;
  n = length(chrs);
  # De-Parallelize
  if(!parallel){foreach::registerDoSEQ()};
  # Apply chromosomal transformations sequentially
  L = foreach(i=1:n) %dopar% {
    # Subset
    key = chrs[i];
    A = X@Data[[key]];
    # Apply
    B = f(A,...);
  }
  names(L) = chrs;
  if(exp.out){
    # Update chromosomal level data
    for(i in 1:n){
      # Subset
      key = chrs[i];
      # Replace
      X@Data[[key]] = L[[key]];
    }
  } else {
    X = L
  }
  # Output
  return(X);
};

#' Apply to FS Experiment
#' 
#' @param X A \code{FSexp} object.
#' @param f Function to apply to each chromosome in \code{FSexp@Data}. Should 
#'   accept objects of class \code{FSchr}.
#' @param exp.out Is the output an object of class \code{FSexp}? Default is 
#'   TRUE.
#' @param parallel Run in parallel? Must register parallel backend first. 
#'   Default is FALSE.
#' @param ... Additional arguments for f.
#' @return List resulting from application of function f to the chromosomes in
#'   \code{FSexp@Data}. If \code{exp.out=T}, the list is formatted as a
#'   \code{FSchr}.
#' 
#' @importFrom foreach "%dopar%" foreach registerDoSEQ
#' @export

apply.FSexp = function(X,f,exp.out=T,parallel=F,...){
  # Chromosomes
  chrs = X@chrs;
  n = length(chrs);
  # De-Parallelize
  if(!parallel){foreach::registerDoSEQ()};
  # Apply chromosomal transformations sequentially
  L = foreach(i=1:n) %do% {
    # Subset
    key = chrs[i];
    A = X@Data[[key]];
    # Apply
    B = f(A,...);
  }
  names(L) = chrs;
  if(exp.out){
    # Update chromosomal level data
    for(i in 1:n){
      # Subset
      key = chrs[i];
      # Replace
      X@Data[[key]] = L[[key]];
    }
  } else {
    X = L
  }
  # Output
  return(X);
};

########################
# Subset Methods
########################

#' subset.DCSchr
#' 
#' Allows subsetting of a DCSchr either using a logical vector, or a maximal
#' allowed pixel separation.
#' 
#' @param x A \code{DCSchr} object.
#' @param ... Unused.
#' @param pix Retain pixels indicated, logical vector.
#' @param dmax Retain pixels with separation not exceeding dmax, numeric.
#' 
#' @export
#' 
#' @examples
#' # Subset pixels on chromosome 1 with maximal separation 10
#' subset(x=DCS::Z@Data[["1"]],dmax=10);

subset.DCSchr = function(x,...,pix,dmax){
  A = x;
  if(!xor(missing(pix),missing(dmax))){stop("Provide either a logical vector, or a maximal pixel separation, not both.")}
  if(!missing(pix)){
    if(!is.logical(pix)){stop("Pix must have class logical.")};
    if(length(pix)!=nrow(A@Coord)){stop("Pix must have length == nrow(A@Coord).")};
    A@Coord = A@Coord[pix,];
    rownames(A@Coord) = NULL;
    A@Counts = A@Counts[pix,];
    rownames(A@Counts) = NULL;
    A@Stats = A@Stats[pix,];
    rownames(A@Stats) = NULL;
  } 
  if(!missing(dmax)){
    if(!is.numeric(dmax)){stop("dmax must have class numeric.")};
    pix = (A@Coord[,"d"]<=dmax);
    A@Coord = A@Coord[pix,];
    rownames(A@Coord) = NULL;
    A@Counts = A@Counts[pix,];
    rownames(A@Counts) = NULL;
    A@Stats = A@Stats[pix,,drop=F];
    rownames(A@Stats) = NULL;
  }
  return(A);
};

#' subset.DCSexp
#' 
#' Allows subsetting of a DCSexp by chromosome and/or maximal pixel separation.
#' 
#' @param x A \code{DCSexp} object. 
#' @param ... Unused.
#' @param chrs Chromosomes to retain, character.
#' @param dmax Maximal pixel separation, numeric.
#' 
#' @export 
#' 
#' @examples 
#' # Subset pixels on chromosome 1 with maximal separation 10
#' subset(x=DCS::Z,chrs="1",dmax=10);

subset.DCSexp = function(x,...,chrs,dmax){
  X = x;
  if(missing(chrs)){chrs = X@chrs};
  # Subset chromosomes 
  X@Data = subset(X@Data,names(X@Data) %in% chrs);
  X@chrs = names(X@Data);
  if(!missing(dmax)){
    # Subset by distance
    if(!is.numeric(dmax)){stop("dmax must have class numeric.")};
    X = apply.DCSexp(X=X,f=subset.DCSchr,dmax=dmax);
  }
  return(X);
};

#' subset.FSchr
#' 
#' Allows subsetting an FSchr using a logical vector and/or the group assignment.
#' 
#' @param x An \code{FSchr} object.
#' @param ... Unused
#' @param foci Foci to reain, logical vector.
#' @param groups Groups to retain, character or factor.
#' 
#' @export 
subset.FSchr = function(x,...,foci,groups){
  B = x;
  if(!missing(foci)){
    if(!is.logical(foci)){stop("Foci must have class logical.")};
    if(length(foci)!=length(B@foci)){stop("Foci must have length == length(B@foci).")};
    B@foci = B@foci[foci];
    B@Stats = plyr::llply(.data=B@Stats,.fun=function(X){X[foci,,drop=F]});
  }
  if(!missing(groups)){
    if(!(is.factor(groups)|is.character(groups))){stop("Groups must have class character or factor.")}
    B@Stats = B@Stats[names(B@Stats) %in% groups];
    B@groups = factor(groups);
  }
  return(B);
};

#' subset.FSexp
#' 
#' Allows subsetting an FSexp by chromosome and/or group.
#' 
#' @param x FSexp.
#' @param ... Unused.
#' @param chrs Chromosomes to retain, character.
#' @param groups Groups to retain, character or factor.
#' 
#' @export 
subset.FSexp = function(x,...,chrs,groups){
  X = x;
  if(!missing(chrs)){
    # Subset chromosomes 
    X@Data = subset(X@Data,names(X@Data) %in% chrs);
    X@chrs = chrs;
  }
  if(!missing(groups)){
    X@Data = plyr::llply(.data=X@Data,.fun=subset.FSchr,groups=groups);
    X@groups = factor(groups);
  }
  return(X);
};

########################
# Print Methods
########################

#' print.DCSchr
#' 
#' @param x A \code{DCSchr} object.
#' @param ... Unused.
#' @param n Number of observations to print. 
#' 
#' @importFrom utils head
#' @export 

print.DCSchr = function(x,...,n=6){
  A = x;
  cat(paste0("Chromosome:  \t ",A@chr,"\n"));
  cat(paste0(c("Samples: \t",paste(A@samples,collapse = "  "),"\n")));
  cat(paste0(c("Groups: \t",paste(as.character(A@groups),collapse = "  "),"\n")));
  cat(paste0(c("Pixels: \t",nrow(A@Coord),"\n\n")));
  cat("Data: \n")
  B = cbind(head(A@Coord,n=n),head(A@Counts,n=n));
  if(ncol(A@Stats)>0){B = cbind(B,head(A@Stats,n=n))};
  print(B);
  cat("\n");
}

#' print.DCSexp
#' 
#' @param x A \code{DCSexp} object.
#' @param ... Unused.
#' @param n Number of observations to print.
#' @param chrs Chromosomes to print, character. 
#' 
#' @export

print.DCSexp = function(x,...,n=6,chrs){
  X = x;
  if(!missing(chrs)){
    if(!is.character(chrs)){stop("Character input expected for chrs.")};
    S = subset.DCSexp(x=X,chrs=chrs);
  } else {
    S = X;
  }
  for(i in 1:length(S@chrs)){
    print.DCSchr(x=S@Data[[i]],n=n);
  }
};

#' print.FSchr
#' 
#' @param x An \code{FSchr} object.
#' @param ... Unused.
#' @param n Number of observations to print.
#' 
#' @importFrom plyr llply
#' @export 

print.FSchr = function(x,...,n=6){
  B = x;
  cat(paste0("Chromosome:  \t ",B@chr,"\n"));
  cat(paste0(c("Groups: \t",paste(as.character(B@groups),collapse = "  "),"\n")));
  cat(paste0(c("Foci: \t",length(B@foci),"\n\n")));
  cat("Stats: \n");
  # Number per group
  npg = ceiling(n/length(B@groups));
  C = do.call(rbind,plyr::llply(.data=B@Stats,.fun=head,n=npg));
  C = cbind("Group"=rep(B@groups,each=npg),"Focus"=rep(B@foci[1:npg],times=length(B@groups)),C);
  rownames(C) = NULL;
  print(C);
  cat("\n");
};

#' print.FSexp
#' 
#' @param x An \code{FSexp} object.
#' @param ... Unused.
#' @param chrs Chromosomes to print, character. 
#' @param n Number of observations to print.
#' 
#' @export
print.FSexp = function(x,...,chrs,n=6){
  X = x;
  if(!missing(chrs)){
    if(!is.character(chrs)){stop("Character input expected for chrs.")};
    S = subset.FSexp(x=X,chrs=chrs);
  } else {
    S = X;
  };
  for(i in 1:length(S@chrs)){
    print.FSchr(x=S@Data[[i]],n=n);
  };
};

#######################
# Merge Method
########################

#' Extract Data from DCSchr
#' 
#' Extracts the coordinates, counts, and stats from a DCS chromosome
#' @param M A \code{DCSchr} object.

extractData = function(M){
  if(class(M)!="DCSchr"){stop("Expecting object of class DCSchr.")};
  Out = data.frame(M@Coord,M@Counts);
  colnames(Out)[4:(4+ncol(M@Counts)-1)] = colnames(M@Counts);
  if(nrow(M@Stats)==nrow(M@Coord)){Out = data.frame(Out,M@Stats)};
  return(Out);
}

#' merge.DCSchr
#' 
#' @param x A \code{DCSchr}. 
#' @param y Another \code{DCSchr}.
#' @param ... Unused.
#' 
#' @return A \code{DCSchr}.
#' 
#' @importFrom methods new
#' @export 

merge.DCSchr = function(x,y,...){
  # Check chromosomes match
  chrs.equal = all.equal(x@chr,y@chr);
  if(chrs.equal!=T){stop("Chromosomes do not match.")}
  chrs = sort(x@chr);
  # Sample sizes
  n1 = length(x@samples);
  n2 = length(y@samples);
  # Sample names
  smps = c(x@samples,y@samples);
  if(length(unique(smps))!=(n1+n2)){stop("Unique sample names required.")};
  # Groups
  groups = as.factor(c(as.character(x@groups),as.character(y@groups)));
  if(length(groups)==length(smps)){groups = groups[order(smps)]};
  smps = sort(smps);
  # Extract data
  A = extractData(x);
  B = extractData(y);  
  # Merge
  C = merge(x=A,y=B,all=T,by=c("i","j","d"));
  C = C[order(C$i,C$j),];
  # Output
  Coord = data.matrix(C[,c("i","j","d")]);
  Counts = data.matrix(C[,smps]);
  nc = ncol(Coord) + ncol(Counts);
  if(ncol(C)>(nc)){Stats = C[,c((nc+1):ncol(C)),drop=F]} else {Stats = data.frame()};
  Out = new("DCSchr",chr=chrs,Coord=Coord,Counts=Counts,
            groups=groups,samples=smps,Stats=Stats);
  return(Out);
}

########################
# Show Methods
########################

#' show for DCSchr
#' @param object A \code{DCSchr} object.
#' @rdname DCSchr-method
#' @importFrom methods show
setMethod(f="show",signature=c(object="DCSchr"),definition=function(object){print.DCSchr(x=object)});

#' show for DCSexp
#' @param object A \code{DCSexp} object.
#' @rdname DCSexp-method
#' @importFrom methods show
setMethod(f="show",signature=c(object="DCSexp"),definition=function(object){print.DCSexp(x=object)});

#' show for FSchr
#' @param object An \code{FSchr} object.
#' @rdname FSchr-method
#' @importFrom methods show
setMethod(f="show",signature=c(object="FSchr"),definition=function(object){print.FSchr(x=object)});

#' show for FSexp
#' @param object An \code{FSexp} object.
#' @rdname FSexp-method
#' @importFrom methods show
setMethod(f="show",signature=c(object="FSexp"),definition=function(object){print.FSexp(x=object)});

########################
# Head Methods
########################

#' head.DCSchr
#' 
#' @param x A \code{DCSchr} object.
#' @param ... Unused.
#' @param n Number of observations to print.
#' 
#' @export 

head.DCSchr = function(x,...,n=6){print.DCSchr(x=x,n=n)};

#' head.DCSexp
#' 
#' @param x A \code{DCSexp} object.
#' @param ... Unused.
#' @param n Number of observations to print.
#' 
#' @export 

head.DCSexp = function(x,...,n=6){print.DCSexp(x=x,n=n)};

#' head.FSchr
#' 
#' @param x An \code{FSchr} object.
#' @param ... Unused.
#' @param n Number of observations to print.
#' 
#' @export 

head.FSchr = function(x,...,n=6){print.FSchr(x=x,n=n)};

#' head.FSexp
#' 
#' @param x An \code{FSexp} object.
#' @param ... Unused.
#' @param n Number of observations to print.
#' 
#' @export 

head.FSexp = function(x,...,n=6){print.FSexp(x=x,n=n)};

########################
# Dim Method
########################

#' dim.DCSchr
#' 
#' @param x A \code{DCSchr} object.
#' 
#' @return Numeric vector, with (number of pixels, number of samples).
#' 
#' @export

dim.DCSchr = function(x){
  A = x;
  c(nrow(A@Coord),length(A@samples))
};

#' dim.DCSexp
#' 
#' @param x A \code{DCSexp} object.
#' 
#' @return A data.frame, with (chromosome, pixels, max pixel separation).
#' 
#' @importFrom plyr llply
#' @export 

dim.DCSexp = function(x){
  X = x;
  A = do.call(rbind,plyr::llply(.data=X@Data,.fun=dim.DCSchr));
  # dmax
  aux = function(A){return(max(A@Coord[,"d"]))};
  d = do.call(rbind,plyr::llply(.data=X@Data,.fun=aux));
  B = data.frame(X@chrs,A[,1],d);
  colnames(B) = c("Chr","Pixels","dMax");
  rownames(B) = NULL;
  return(B);
};

#' dim.FSchr
#' 
#' @param x An \code{FSchr} object.
#' 
#' @return A numeric vector, with (number of foci).
#' 
#' @export

dim.FSchr = function(x){
  B = x;
  length(B@foci)
};

#' dim.FSexp
#' 
#' @param x A \code{FSexp} object.
#' 
#' @return A data.frame, with (chromosomes, number of foci).
#' 
#' @importFrom plyr llply
#' @export
 
dim.FSexp = function(x){
  X = x;
  A = do.call(rbind,plyr::llply(.data=X@Data,.fun=dim.FSchr));
  # dmax
  B = data.frame("Chr"=X@chrs,"Foci"=A[,1]);
  rownames(B) = NULL;
  return(B);
};

########################
# as.data.frame Methods
########################

#' as.data.frame.DCSchr
#' 
#' @param x A \code{DCSchr} object.
#' @param row.names Unused.
#' @param optional Unused.
#' @param ... Unused.
#' 
#' @export

as.data.frame.DCSchr = function(x,row.names,optional,...){
  A = x;
  Out = data.frame("Chr"=A@chr,A@Coord,A@Counts,A@Stats);
  rownames(Out) = NULL;
  return(Out);
};

#' as.data.frame.DCSexp
#' 
#' @param x A \code{DCSexp} object.
#' @param row.names Unused.
#' @param optional Unused.
#' @param ... Unused.
#'  
#' @importFrom foreach "%do%" foreach
#' @export 

as.data.frame.DCSexp = function(x,row.names,optional,...){
  X = x;
  i = NULL;
  Out = foreach(i=1:length(X@chrs),.combine=rbind) %do% {
    return(as.data.frame.DCSchr(X@Data[[i]]));
  }
  row.names(Out) = NULL;
  return(Out);
};

#' as.data.frame.FSchr
#' 
#' @param x An \code{FSchr} object.
#' @param row.names Unused.
#' @param optional Unused.
#' @param ... Unused.
#' 
#' @importFrom foreach "%do%" foreach
#' @export

as.data.frame.FSchr = function(x,row.names,optional,...){
  B = x;
  i = NULL;
  Out = foreach(i=1:length(B@groups),.combine=rbind) %do% {
    return(data.frame("Chr"=B@chr,"Group"=B@groups[i],"Focus"=B@foci,B@Stats[[i]]));
  }
  rownames(Out) = NULL;
  return(Out);
};

#' as.data.frame.FSexp
#' 
#' @param x An \code{FSexp} object. 
#' @param row.names Unused.
#' @param optional Unused.
#' @param ... Unused.
#'  
#' @importFrom foreach "%do%" foreach
#' @export

as.data.frame.FSexp = function(x,row.names,optional,...){
  X = x;
  i = NULL;
  Out = foreach(i=1:length(X@chrs),.combine=rbind) %do% {
    return(as.data.frame.FSchr(X@Data[[i]]));
  }
  row.names(Out) = NULL;
  return(Out);
};

########################
# Summary Method
########################

#' summary.DCSexp
#' 
#' @param object A \code{DCSexp} object.
#' @param ... Unused.
#' @return A list containing 1. the dimensions of each chromosome in the 
#'   experiment, and 2. the group assignments.
#'   
#' @export

summary.DCSexp = function(object,...){
  X = object;
  A = dim.DCSexp(X);
  print(A);
  cat("\n")
  B = data.frame("Sample"=X@samples,"Group"=X@groups);
  print(B);
  return(list("Dimensions"=A,"Groups"=B));
};