# Purpose : Using the border score, identify candidate boundaries in each group.

#' Call Peaks
#' 
#' Scans across a chromosome. Calculates the mean and variance of the border
#' score within the window [focus-h,focus+h]. Identifies foci with a
#' standardized border score above the threshold.
#' @param foci Focus coordinate.
#' @param b Border score.
#' @param h Window half-width.
#' @param t Border score threshold.
#' @importFrom plyr llply
#' @importFrom stats var

callPeaks = function(foci,b,h,t){
  # Identify candidate peaks
  aux = function(y){
    # Positions to include
    keep = (abs(foci-y) <= h)
    # Local mean and variance
    mu = mean(b[keep]);
    v = var(b[keep]);
    # Z scores
    z = (b[keep] - mu)/sqrt(v);
    # Out
    Out = (foci[keep])[z>=t];
  }
  peaks = plyr::llply(.data=foci,.fun=aux);
  peaks = sort(Reduce(f=union,x=peaks));
  return(peaks);
}

#' Thin Peaks
#' 
#' If multiple candidate peaks are within the minimum separation distance,
#' evaluates the border score at each, and retains the representative peak with
#' the highest border score.
#' @param p Peak coordinates
#' @param b Peak border scores
#' @param w Minimum separation
#' @importFrom plyr aaply

thinPeaks = function(p,b,w){
  flag = T;
  peaks = p;
  b.p = b;
  while(flag){
    # Reduce to maximal peaks
    aux = function(y){
      # Positions to include 
      keep = (abs(y-peaks)<=w);
      # Maximum
      topPeak = which.max(b.p[keep]);
      topPeak = (peaks[keep])[topPeak];
      return(topPeak);
    }
    peaks = plyr::aaply(.data=peaks,.fun=aux,.margins=1);
    peaks = unique(peaks);
    b.p = b[p %in% peaks];
    if(length(peaks)<2){break};
    flag = (min(abs(diff(peaks)))<=w);
  }
  # Output
  return(peaks);
}

#' Candidate Differential Boundaries
#' 
#' Scans across each chromosome. Calculates the local mean and variance of the border
#' score within the window [focus-h,focus+h]. Identifies foci with a border score at least
#' \emph{t} standard deviations above the mean. If multiple foci are within distance w of 
#' one another, the representative focus with the greatest border score is retained. 
#' @param B FS experiment produced by \code{BorderScores}.
#' @param h Window half-width. 
#' @param w Minimum boundary separation.
#' @param t Threshold local Z-score for candidate boundaries.
#' @param parallel Run in parallel? Must register parallel backend first. 
#' @return data.frame of candidate boundaries
#' @importFrom foreach "%do%" "%:%" foreach registerDoSEQ
#' @export

callBoundaries = function(B,h,w=10,t=2,parallel=F){
  # De-Parallelize
  if(!parallel){foreach::registerDoSEQ()};
  # Groups
  ng = length(B@groups);
  # Chromosomes
  nc = length(B@chrs);
  # Loop over groups
  i = j = NULL;
  Out = foreach(i=1:ng,.combine=rbind) %:% 
    foreach(j=1:nc,.combine=rbind) %dopar% {
      # Subset chromosome
      C = B@Data[[B@chrs[j]]];
      # Foci
      foci = C@foci;
      # Border scores
      Stats = C@Stats[[B@groups[i]]];
      border = Stats$Border;
      # Remove NAs
      keep = !is.na(border);
      Stats = Stats[keep,];
      border = border[keep];
      foci = foci[keep];
      # Candidate peaks
      Peaks = callPeaks(foci=foci,b=border,h=h,t=t);
      if(length(Peaks)>0){
        Peaks = thinPeaks(p=Peaks,b=border[foci %in% Peaks],w=w);
        Peaks = data.frame("Chr"=B@chrs[j],"Group"=B@groups[i],
                           "Focus"=Peaks,Stats[foci %in% Peaks,]); 
      } else {Peaks = NULL};
      # Output
      return(Peaks);
    }
  rownames(Out) = NULL;
  return(Out);
}