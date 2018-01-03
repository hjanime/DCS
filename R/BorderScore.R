#' @include ClassDefinitions.R
NULL

# Purpose : Functions to calculate a border score for each focus in the genome.
# Border scores compare the signal intensity in a candidate inter-TAD region with
# the mean border score in the flanking regions. 

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
#' @importFrom plyr llply

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
#' @importFrom plyr llply

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
#' @importFrom plyr llply

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
#' @param G Chromosome level data structured as (i,j,groupMean)
#' @param agg FUNCTION for signal aggregation within each region. Defaults to sum.
#' @return A matrix containing the focus, region summary statistics, and the border
#' score. Regions (A,B,C) refer to the upstream, downstream, and inter-TAD regions 
#' respectively. Signal is the aggregated signal in the region. Pixels is the total 
#' number of pixels in the region, and prop is the proportion of pixels which 
#' were populated. 
#' @importFrom plyr alply

chrBorderScore = function(G,agg=sum){
  
  # Restrict to positions where i = j
  Foci = G[G[,"i"]==G[,"j"],];
  # Upper limit
  l = min(G[,"i"]);
  u = max(G[,"j"]);
  # Maximum Separation
  dmax = max(G[,"d"])
  # Function to calculate border score; r is a row of G.
  aux = function(r){
    # Focus
    current.focus = as.numeric(r[1]);
    # Null vector, in case one of the signal regions is empty
    nullvec = c(current.focus,rep(NA,times=10));
    # Upstream TAD
    A = triCoordLeft(f=current.focus,d=dmax,l=l,u=u);
    mA = nrow(A);
    if(is.null(mA)){return(nullvec)};
    # Downstream TAD
    B = triCoordRight(f=current.focus,d=dmax,l=l,u=u);
    mB = nrow(B);
    if(is.null(mB)){return(nullvec)};
    # Inter-TAD Region
    C = triCoord(f=current.focus,d=dmax,l=l,u=u);
    mC = nrow(C);
    if(is.null(mC)){return(nullvec)};
    # Merge in statistic
    colnames(A) = colnames(B) = colnames(C) = c("i","j");
    # Merge in statistic
    A = merge(x=A,y=G,by=c("i","j"));
    B = merge(x=B,y=G,by=c("i","j"));
    C = merge(x=C,y=G,by=c("i","j"));
    # Occupied cells
    oA = nrow(A);
    oB = nrow(B);
    oC = nrow(C);
    # Aggregate
    if(nrow(A)==0){sA=0} else {sA = (1/mA)*agg(A[,"groupMean"]);};
    if(nrow(B)==0){sB=0} else {sB = (1/mB)*agg(B[,"groupMean"]);};
    if(nrow(C)==0){sC=0} else {sC = (1/mC)*agg(C[,"groupMean"]);};
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
#' @param A Chromosome-level data from HiC experiment.
#' @param agg FUNCTION for signal aggregation within each region. Defaults to sum.
#' @importFrom foreach "%do%" foreach

borderWrap = function(A,agg=sum){
  # Samples
  s = A@samples;
  ns = length(s);
  # Group levels
  g = A@groups;
  L = sort(unique(g));
  ng = length(L);
  # Foci
  i = j = NULL; 
  foci = A@Coord[A@Coord[,"i"]==A@Coord[,"j"],"i"];
  # Loop over groups
  B = foreach(i=1:ng) %do% {
    # Counts
    Counts = A@Counts;
    # Subset group
    Sub = Counts[,g==L[i],drop=F];
    if(ncol(Sub)>1){Sub = apply(Sub,MARGIN=1,FUN=mean)};
    # Bind group mean to coordinates
    G = cbind(A@Coord,Sub);
    colnames(G)[4] = "groupMean";
    # Calculate border scores
    C = data.frame(chrBorderScore(G=G,agg=agg));
    # Quality Metrics
    C$minPix = apply(C[,c("pixelsA","pixelsB","pixelsC")],MARGIN=1,FUN=min);
    C$minProp = apply(C[,c("propA","propB","propC")],MARGIN=1,FUN=min);
    C = C[,c("Border","minPix","minProp")];
    return(C);
  }
  names(B) = L;
  # Format as FSchr
  Out = new(Class="FSchr",chr=A@chr,foci=foci,groups=L,Stats=B)
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
#' @param agg Function used to aggregate signal within each region. 
#' @param parallel Run in parallel? Must register parallel backend first. 
#' @return Returns a list with one data.frame per group.
#' @importFrom methods new
#' @export

BorderScores = function(X,agg=sum,parallel=F){
  # Obtain border scores for each chromosome
  B = apply.DCSexp(X=X,f=borderWrap,exp.out=F,parallel=parallel,agg=agg);
  # Bind into FSexp
  Out = new(Class="FSexp",chrs=X@chrs,groups=B[[1]]@groups,Data=B);
  return(Out);
}