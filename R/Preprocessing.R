########################
# 1.Impute Missing Pixels to Zero 
########################

#' Coordinates of Matrix Diagonals
#'
#' @param a Stat of sequence
#' @param b End of sequence
#' @param d Maximum pixel separation
#' @importFrom plyr llply

genDiag = function(a,b,d){
  # Separations
  sep = seq(from=0,to=d);
  # Coords
  aux = function(x){
    i = seq(from=a,to=(b-x));
    j = seq(from=(a+x),to=b);
    return(cbind(i,j));
  }
  Out = plyr::llply(.data=sep,.fun=aux);
  Out = do.call(rbind,Out);
  Out = cbind(Out,Out[,2]-Out[,1]);
  return(Out);
}

#' Coordinates of Missing Pixels
#' 
#' Enumerates coordinates of pixels NOT present in the data.
#' @param Coords Matrix of pixel coordinates

missedCoord = function(Coords){
  # Possible cells
  allC = genDiag(a=min(Coords[,1]),b=max(Coords[,2]),d=max(Coords[,3]));
  # Finds the coordinates of pixels in allC and not in C
  Out = rbind(allC,Coords);
  Out = Out[!duplicated(Out, fromLast=TRUE) & seq(nrow(Out)) <= nrow(allC),];
  Out = matrix(Out,ncol=3);
  return(Out);
}

#' Impute Missing Pixels
#' 
#' @param Counts Matrix of pixel intensities
#' @param Coords Matrix of pixel coordinates

impute = function(Counts,Coords){
  # Missing coords
  Miss = missedCoord(Coords);
  if(nrow(Miss)==0){
    Out = list("Impute"=F);
  } else {
    # Impute counts to zero
    Zero = array(0,dim=c(nrow(Miss),ncol(Counts)));
    # Imputed Counts
    Counts.out = rbind(Counts,Zero);
    rownames(Counts.out) = NULL;
    # Imputed Coords
    Coords.out = rbind(Coords,Miss);
    rownames(Coords.out) = NULL;
    # Output
    Out = list("Counts"=Counts.out,"Coords"=Coords.out,"Impute"=T);
  }
  return(Out);
}

#' Wrapper for Imputing Missing Pixels
#' 
#' @param A DCSchr

imputeWrap = function(A){
  # Counts
  Counts = A@Counts;
  # Impute NA to zero
  Counts[is.na(Counts)] = 0;
  A@Counts = Counts;
  # Coords
  Coords = A@Coord;
  # Impute
  imputeOut = impute(Counts,Coords);
  if(imputeOut$Impute==F){
    # If no pixels were imputed, return the original DCSchr
    return(A);
  } else {
    A@Counts = imputeOut$Counts;
    A@Coord = imputeOut$Coords;
    return(A);
  }
}

########################
# 2. Remove Empty Foci 
########################

#' Remove Pixels Corresponding to Empty Foci
#' 
#' @param A Chromosome level data
#' @importFrom plyr aaply 

pruneEmpty = function(A){
  # Combine Counts and Coords
  B = cbind(A@Coord,A@Counts);
  # Split data into foci and off-diagonal pixels
  B = B[order(B[,"d"]),];
  Foci = B[B[,"d"]==0,];
  # Coords and Counts
  Coord = Foci[,c("i","j")];
  Counts = Foci[,c(4:ncol(B))]
  # Positions to remove
  aux = function(x){sum(as.numeric(x)==0)>0};
  Remove = plyr::aaply(.data=Counts,.margins=1,.fun=aux);
  # Foci to remove
  Coord.r = Coord[Remove,1];
  # Check: Count.r = Count.f[Remove,];
  # Positions to remove
  key = (B[,"i"] %in% Coord.r)|(B[,"j"] %in% Coord.r);
  # Positions retained
  Retained = B[!key,];
  rownames(Retained) = NULL;
  # Output
  A@Coord = Retained[,c(1:3)];
  A@Counts = Retained[,c(4:ncol(Retained))];
  return(A);
}

########################
# 3. Pre-processing function
########################

#' Wrapper Function for Pre-processing
#' 
#' @param A Chromosome level data
#' @param miss Impute missing pixels to zeros?
#' @param zero Remove pixels corresponding to empty foci? 

wrapPP = function(A,miss=T,zero=T){
  # Output
  Out = A;
  # Impute
  if(miss){Out = imputeWrap(A=Out)}
  # Filter
  if(zero){Out = pruneEmpty(A=Out)};
  return(Out);
}

#' Pre-process HiC Experiment
#'
#' @param X HiC experiment.
#' @param miss Impute missing pixels to zeros?
#' @param zero Remove pixels corresponding to empty foci? 
#' @param parallel Run in parallel? Must register parallel backend first.  
#' @export

preProcess = function(X,miss=T,zero=T,parallel=F){
  Out = apply.DCSexp(X=X,f=wrapPP,miss=miss,zero=zero,parallel=parallel);
  return(Out);
}

########################
# 4. Log fold change filter
########################

#' Calculate log2 Fold Change
#' 
#' @param Counts Pixel intensities
#' @param g Sample group assignments
#' @importFrom foreach "%do%" foreach

calcLFC = function(Counts,g){
  # Group levels
  L = sort(unique(g));
  # Number of groups
  ng = length(L);
  # Group means
  i = NULL;
  Means = foreach(i=1:ng,.combine=cbind) %do% {
    # Subset group columns
    Sub = Counts[,g %in% L[i],drop=F]
    # Collapse
    if(ncol(Sub)>1){
      Sub = apply(Sub,MARGIN=1,FUN=mean,na.rm=T)
    }
    return(Sub);
  }
  # Reference group
  Ref = Means[,1];
  # Group log fold changes
  j = NULL;
  LFC = foreach(j=2:ng,.combine=cbind) %do% {
    # Ratio
    Rho = (Means[,j]+1)/(Ref+1);
    # Calculate log fold change
    Out = log2(Rho);
    return(Out)
  }
  Out = as.matrix(LFC);
  # Output
  colnames(Out) = paste0("lfc_",L[2:ng],"v",L[1]);
  return(Out);
}

#' Filter on Log Fold Change
#' 
#' @param A Chromosome level data
#' @param aT Absolute threshold 
#' @param qT Quantile threshold
#' @importFrom stats quantile
#' @importFrom foreach "%do%" foreach

# Function to filter based on log fold change
wrapLfcFilter = function(A,aT,qT){
  # Combine Counts and Coords
  B = cbind(A@Coord,A@Counts);
  # Split data into foci and off-diagonal pixels
  B = B[order(B[,"d"]),];
  Foci = B[B[,"d"]==0,];
  # Log fold changes
  LFC = calcLFC(Counts=Foci[,4:ncol(Foci)],g=A@groups);
  # Thresholding
  if(!missing(qT)){
    # Quantile thresholding
    aux = function(x){quantile(x,probs=c(qT,1-qT))};
    tau = apply(LFC,MARGIN=2,FUN=aux);
  } else {
    # Absolute thresholding
    tau = rep(aT,ncol(LFC));
    tau = rbind(-tau,tau);
  }
  # Foci to keep
  i = j = NULL;
  Remove = foreach(i=1:ncol(tau),.combine=cbind) %do% {
    # Subset 
    tau.sub = tau[,i];
    LFC.sub = LFC[,i];
    # Removal indicator
    Out = (LFC.sub < tau.sub[1])|(LFC.sub > tau.sub[2]);
    return(Out);
  }
  Remove = as.matrix(Remove);
  # Reduce
  Remove = (apply(Remove,MARGIN=1,FUN=sum)>0);
  # Foci to discard
  Coord.r = Foci[Remove,"i"];
  # Positions to remove
  key = (B[,"i"] %in% Coord.r)|(B[,"j"] %in% Coord.r);
  # Retained
  Retained = B[!key,];
  # Output
  A@Coord = Retained[,c(1:3)];
  A@Counts = Retained[,c(4:ncol(Retained))];
  return(A);
}

########################
# 5. Add LFC to an experiment
########################

#' Filter Based on Log Fold Change Relative to Reference Group
#'
#' @param X A \code{DCSexp} object. 
#' @param aT Absolute LFC threshold
#' @param qT Quantile for LFC thresholding
#' @param parallel Run in parallel? Must register parallel backend first. 
#' @export

filterLFC = function(X,aT,qT,parallel=F){
  if(!xor(missing(aT),missing(qT))){stop("Provide either an absolute threshold, or a quantile on which to threshold.")};
  Out = apply.DCSexp(X=X,f=wrapLfcFilter,aT=aT,qT=qT,parallel=parallel);
  return(Out);
}

addLfcWrap = function(A){
  # Calculate LFC
  LFC = calcLFC(Counts=A@Counts,g=A@groups);
  rownames(LFC) = NULL;
  # Store LFC stat
  l = nrow(A@Stats);
  if(l==0){A@Stats = data.frame("lfc"=LFC)
    } else if(l>0){A@Stats = cbind(A@Stats,"lfc"=LFC)}; 
  # Output
  return(A);
}

#' Add Log Fold Change Relative to Reference Group
#'
#' Note: Add LFC after normalization. 
#' @param X A \code{DCSexp} object.
#' @param parallel Run in parallel? Must register parallel backend first. 
#' @export

addLFC = function(X,parallel=F){
  Out = apply.DCSexp(X=X,f=addLfcWrap,parallel=parallel);
  return(Out);
}