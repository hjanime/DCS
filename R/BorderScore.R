#' Purpose : Functions to calculate a border score for each focus in the genome.
#' Border scores compare the signal intensity in a putative inter-TAD region with
#' the mean border score in the flanking regions. 

#' Enumerate Inter-Tad Region
#' 
#' Enumerates pixels in an upper triangular window above a focus. Maximal 
#' pixel separation d corresond to the base-width of the triangle in pixels. 
#' Assuming the pixels are indexed (i,j) with i <= j, the lower bound l = min(i), 
#' while the upper bound u = max(j). 
#' @param f Focus. 
#' @param d Maximum pixel separation. 
#' @param l Lower bound on pixel indices. 
#' @param u Upper bound on pixel indices.  
#' @import plyr

triCoord = function(f,d,l,u){
  # f : Focus
  # d : Max bin separation
  # u : Upper limit
  
  # Function to enumerate coordinates (down columns);
  aux = function(z){
    if(f+z<=u){
      a = seq(from=max(l,f-d+z),to=f);
      b = rep(x=f+z,times=length(a));
      A = cbind(a,b);
      return(A);
    }
  };
  Out = plyr::llply(.data=seq(from=0,to=d),.fun=aux);
  Out = do.call(rbind,Out);
  # Output
  return(Out);
}

#' Enumerate Upstream Intra-TAD Region
#' 
#' Enumerates pixels in an upper triangular window upstream of a focus. Maximal 
#' pixel separation d corresond to the base-width of the triangle in pixels. 
#' Assuming the pixels are indexed (i,j) with i <= j, the lower bound l = min(i), 
#' while the upper bound u = max(j). 
#' @param f Focus. 
#' @param d Maximum pixel separation. 
#' @param l Lower bound on pixel indices. 
#' @param u Upper bound on pixel indices. 
#' @import plyr 

triCoordLeft = function(f,d,l,u){
  # f : Focus
  # r : Resolution
  # d : Max bin separation
  # u : Upper limit
  
  # Upper left coord of region
  UL = max(l,f-d-1);
  # Function to enumerate coordinates (down columns);
  aux = function(z){
    if((UL+z)<f){
      a = seq(from=UL,to=(UL+z));
      b = rep(x=(UL+z),times=length(a));
      A = cbind(a,b);
      return(A);
    }
  }
  Out = plyr::llply(.data=seq(from=0,to=d),.fun=aux);
  Out = do.call(rbind,Out);
  # Output
  return(Out);
}

#' Enumerate Downstream Intra-TAD Region
#' 
#' Enumerates pixels in an upper triangular window downstream of a focus. Maximal 
#' pixel separation d corresond to the base-width of the triangle in pixels. 
#' Assuming the pixels are indexed (i,j) with i <= j, the lower bound l = min(i), 
#' while the upper bound u = max(j). 
#' @param f Focus. 
#' @param d Maximum pixel separation. 
#' @param l Lower bound on pixel indices. 
#' @param u Upper bound on pixel indices.  
#' @import plyr

triCoordRight = function(f,d,l,u){
  # Upper Left of Region
  UL = f+1;
  # Function to enumerate coordinates (down columns);
  aux = function(z){
    if((UL+z)<= u){
      a = seq(from=UL,to=UL+z);
      b = rep(x=(UL+z),times=length(a));
      A = cbind(a,b);
      return(A);
    }
  }
  Out = plyr::llply(.data=seq(from=0,to=d),.fun=aux);
  Out = do.call(rbind,Out);
  # Output
  return(Out);
}

#' Chromosome-level Border Score Calculation
#' 
#' For each focus, aggregates signal in the inter-TAD region, and within the
#' upstream and downstream intra-TAD regions. Calculates a border score as the 
#' ratio of mean signal in the intra-TAD regions to mean signal in the inter-TAD region.
#' @param W Chromosome level matrix structured as (i,j,Stat);
#' @param agg FUNCTION for signal aggregation within each region. Defaults to sum.
#' @return A matrix containing the focus, region summary statistics, and the border
#' score. Regions (A,B,C) refer to the upstream, downstream, and inter-TAD regions 
#' respectively. Signal is the aggregated signal in the region. Pixels is the total 
#' number of pixels in the region, and prop is the proportion of pixels which 
#' were populated. 
#' @import plyr

chrBorderScore = function(W,agg=sum){
  
  # Restrict to positions where i = j
  Foci = W[W[,1]==W[,2],];
  # Upper limit
  l = min(W[,1]);
  u = max(W[,2]);
  # Separation
  d = max(abs(W[,1]-W[,2]));
  # Function to calculate border score
  # x is a row of W
  aux = function(x){
    # Focus
    current.focus = as.numeric(x[1]);
    # Upstream TAD
    A = triCoordLeft(f=current.focus,d=d,l=l,u=u);
    mA = nrow(A);
    if(is.null(mA)){return()};
    # Downstream TAD
    B = triCoordRight(f=current.focus,d=d,l=l,u=u);
    mB = nrow(B);
    if(is.null(mB)){return()};
    # Inter-TAD Region
    C = triCoord(f=current.focus,d=d,l=l,u=u);
    mC = nrow(C);
    if(is.null(mC)){return()};
    # Merge in statistic
    colnames(A) = colnames(B) = colnames(C) = c("i","j");
    # Merge in statistic
    A = merge(x=A,y=W,by=c("i","j"));
    B = merge(x=B,y=W,by=c("i","j"));
    C = merge(x=C,y=W,by=c("i","j"));
    # Occupied cells
    oA = nrow(A);
    oB = nrow(B);
    oC = nrow(C);
    # Aggregate
    if(nrow(A)==0){sA=0} else {sA = (1/mA)*agg(A[,3]);};
    if(nrow(B)==0){sB=0} else {sB = (1/mB)*agg(B[,3]);};
    if(nrow(C)==0){sC=0} else {sC = (1/mC)*agg(C[,3]);};
    # Region A
    outA = c(sA,mA,oA/mA);
    names(outA) = c("signalA","pixelsA","propA");
    # Region B
    outB = c(sB,mB,oB/mB);
    names(outB) = c("signalB","pixelsB","propB");
    # Region C
    outC = c(sC,mC,oC/mC);
    names(outC) = c("signalC","pixelsC","propC");
    # Border
    border = (sA+sB)/(2*sC);
    out = c("Focus"=current.focus,outA,outB,outC,"Border"=border);
    return(out);
  }
  # Output
  Out = plyr::alply(.data=Foci,.margins=1,.fun=aux);
  Out = do.call(rbind,Out);
  return(Out);
};

#' Wrapper Function for Chromosome-level Border Score Calculation
#' 
#' @param W Chromosome-level data from HiC experiment.
#' @param s Sample names.
#' @param g Group assignments. 
#' @param agg FUNCTION for signal aggregation within each region. Defaults to sum.

borderWrap = function(W,s,g,agg){
  # Chromosome
  chr = unique(W$Chr);
  chr2 = unique(W$Chr2);
  # Group levels
  L = sort(unique(g));
  ng = length(L);
  # Loop over groups
  B = foreach(i=1:ng,.combine="rbind") %do% {
    # Samples to subset
    s.sub = s[g==L[i]]
    # Subset
    Sub = W[,s.sub,drop=F];
    if(ncol(Sub)>1){Sub = apply(Sub,MARGIN=1,FUN=mean)};
    Sub = data.frame(W[,c("i","j")],Sub);
    colnames(Sub)[3] = "groupMean";
    # Calculate border scores
    B = data.frame(chrBorderScore(W=Sub,agg=agg));
    B$Group = paste0("g",L[i]);
    return(B);
  }
  Out = B[,c("Focus","Group","Border")];
  # Quality metrics
  Out$minPix = apply(B[,c("pixelsA","pixelsB","pixelsC")],MARGIN=1,FUN=min);
  Out$minProp = apply(B[,c("propA","propB","propC")],MARGIN=1,FUN=min);
  # Format
  Out$Chr = chr;
  Out$Chr2 = chr2;
  Out = Out[,c("Focus","Chr","Chr2","Group","Border","minPix","minProp")]
  return(Out);
};

#' Experiment-level Border Score Calculation
#' 
#' Calculates a border scores by comparing the mean pixel intensity in 
#' triangular windows upstream and downstream of the focus to the mean pixel 
#' intensity in the intervening region.
#' 
#' minPix refers to the minimum number of pixels observed in one of the
#' aggregation regions. minProp refers to the minimum proportion of pixels
#' observed in one of the aggregation regions.
#' 
#' @param X HiC experiment.
#' @param s Sample names.
#' @param g Sample group assignments.
#' @param agg agg FUNCTION for signal aggregation within each region. Defaults
#'   to sum.
#' @param cores Cores to use if running in parallel.
#' @return Returns a list with one data.frame per group.
#' @export

getBorderScores = function(X,s,g,agg=sum,cores=1){
  # Obtain border scores for each chromosome
  A = expApply(X=X,f=borderWrap,cores=cores,s=s,g=g,agg=agg);
  # Collapse experiment
  A = do.call(rbind,A);
  row.names(A) = NULL;
  # Create separate data.frame for each group
  L = sort(unique(A$Group));
  ng = length(L);
  Out = foreach(i=1:ng) %do% {
    return(A[A$Group==L[i],])
  }
  names(Out) = L;
  return(Out);
}