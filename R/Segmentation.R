#' Call Peaks
#' 
#' Scans across a chromosome. Calculates the mean and variance of the
#' border score within the window. Identifies foci with a standardized
#' border score above the threshold. 
#' @param x Focus coordinates.
#' @param b Border score.
#' @param h Window half-width.
#' @param t Border score threshold.
#' @import plyr

callPeaks = function(x,b,h,t){
  # Identify candidate peaks
  aux = function(y){
    # Positions to include
    keep = (abs(x-y) <= h)
    # Local mean and variance
    mu = mean(b[keep]);
    v = var(b[keep]);
    # Z scores
    z = (b[keep] - mu)/sqrt(v);
    # Out
    Out = (x[keep])[z>=t];
  }
  peaks = plyr::llply(.data=x,.fun=aux);
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
#' @import plyr

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
    flag = (min(abs(diff(peaks)))<=w);
  }
  # Output
  return(peaks);
}

#' Segment TADs
#' 
#' For each chromosome, identifies the coordinates of foci that
#' are at least t standard deviations above the local mean. If 
#' multiple foci are within distance w of one another, the
#' representative focus with the greatest border score is retained. 
#' @param chr Chromosome.
#' @param x Forcus coordinates.
#' @param b Border score. 
#' @param h Window half-width. 
#' @param w Minimum boundary separation.
#' @param t Threshold local Z-score for candidate boundaries. 
#' @import doMC
#' @export 

segmentTADs = function(chr,x,b,h,w,t){
  # Data frame
  df = data.frame(chr,x,b);
  # Chromosomes
  chrs = unique(chr);
  nc = length(chrs);
  A = foreach(i=1:nc,.combine=rbind) %do% {
    # Subset chromosome
    Sub = df[df$c==chr[i],2:3];
    # All foci
    all.foci = Sub[,1];
    # All border scores
    all.border = Sub[,2];
    # Candidate peaks
    Peaks = callPeaks(x=all.foci,b=all.border,h=h,t=t);
    if(length(Peaks)>0){
      Peaks = thinPeaks(p=Peaks,b=all.border[all.foci %in% Peaks],w=w);
      Peaks = data.frame("Focus"=Peaks,"Chr"=chrs[i]) 
      } else {Peaks = NULL};
    # Output
    return(Peaks);
  }
  return(A);
}

#' Candidate Differential Boundaries
#' 
#' @param B Matrix containing the following columns:
#' 1 & 2. Chromosome, 3. Focus, 4. Group, 5. Border score
#' @param h Window half-width. 
#' @param w Minimum boundary separation.
#' @param t Threshold local Z-score for candidate boundaries. 
#' @import doMC
#' @export

getBoundaries = function(B,h,w,t){
  # Relabel
  names(B) = c("Chr","Chr2","Focus","Group","Border");
  C = unique(B[,c("Chr","Chr2")]);
  # Groups
  Groups = unique(B$Group);
  ng = length(Groups);
  # Loop over groups
  Out = foreach(i=1:ng) %do% {
    # Subset group
    A = B[B$Group==Groups[i],];
    # Segment TADs
    Seg = segmentTADs(chr=A$Chr,x=A$Focus,b=A$Border,h=h,w=w,t=t);
    # Add chr2
    Seg = merge(x=Seg,y=C,by="Chr");
    # Add group
    Seg = cbind(Seg,"Group"=Groups[i]);
    # Reorder columns
    Seg = Seg[,c("Focus","Chr","Chr2","Group")];
  }
  names(Out) = Groups;
  # Union 
  Union = do.call(rbind,Out);
  Union = Union[order(Union$Chr,Union$Focus),];
  Union = unique(Union);
  rownames(Union) = NULL;
  # Output
  Out$Union = Union;
  return(Out);
}