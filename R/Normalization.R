########################
# 1. Normalization
########################

#' Diagonal Means
#' 
#' Calculate the mean intensity across samples for each pixel separation
#' \emph{d}.
#' @param A DCS chromosome.
#' @importFrom reshape2 melt

diagMeans = function(A){
  # Distance
  d = factor(A@Coord[,"d"]);
  # Levels
  L = sort(unique(d));
  # Intensities
  Counts = data.frame(A@Counts);
  Counts$d = d;
  # Reshape
  R = reshape2::melt(Counts,id.vars="d");
  colnames(R)[2:3] = c("Sample","y");
  # Aggregate
  Out = aggregate(y ~ Sample + d, data=R, FUN=mean);
  Out$Chr = unique(A@chr);
  return(Out);
}

#' Scale Factor Estimation
#' 
#' Function to calculate separation based scale factors. In particular, the scale factor
#' for the \emph{i}th sample is calcuate as the mean across samples at separation \emph{d}, \eqn{m(d)},
#' divided by the mean in sample \emph{i} at separation \emph{d}, \eqn{m_{i}(d)}. 
#' @param X HiC experiment
#' @param parallel Run in parallel? Must register parallel backend first. 
#' 
#' @importFrom foreach "%do%" foreach
#' @importFrom stats aggregate
#' @export 

scaleFactors = function(X,parallel=F){
  # Samples
  s = X@samples;
  # Number of samples
  ns = length(X@samples);
  # Means by separation
  Z = apply.DCSexp(X=X,f=diagMeans,exp.out=F,parallel=parallel);
  Z = do.call(rbind,Z);
  # Calculate sample specific mean at separation d
  y = Sample = d = NULL;
  Z1 = aggregate(y ~ Sample + d, data=Z, FUN=mean);
  # Calculate mean across samples at separation d
  Z2 = aggregate(y ~ d, data=Z1, FUN=mean);
  # Loop over samples
  i = NULL;
  S = foreach(i=1:ns,.combine=cbind) %do% {
    # Subset sample.
    Sub = Z1[Z1$Sample==s[i],];
    # Scale factors (mean across samples / mean in sample i)
    scale = Z2$y/Sub$y;
    # Correct for potential division by zero
    scale[is.nan(scale)] = 0;
    return(scale);
  }
  # Format
  S = data.matrix(S);
  colnames(S) = s;
  S = cbind(S,"d"=as.numeric(as.character(Z2$d)));
  # Output
  return(S);
}

#' Rescale Intensities
#' 
#' Rescales the observed HiC intensities along the diagonals of the contact matrix.
#' @param A DCS chromosome.
#' @param S Matrix of scale factors
diagScale = function(A,S){
  # Samples
  s = A@samples;
  # Number of samples
  ns = length(s);
  # Remove counts
  Counts = cbind(A@Counts,"d"=A@Coord[,"d"]);
  # Merge in scale factors
  B = data.matrix(merge(x=Counts,y=S,by="d"));
  # Separate intensities and scale factors for element-wise multiplication
  B1 = B[,2:(ns+1)]
  B2 = B[,(ns+2):ncol(B)]
  # Scaled Counts
  B3 = B1 * B2;
  colnames(B3) = s;
  # Output
  A@Counts = B3;
  return(A);
}

#' Normalize HiC Samples
#' 
#' Normalizes a HiC experiment s.t. all samples have the same mean pixel intensity at each 
#' possible pixel separation \eqn{d\in\{1,\cdots,d_{\max}\}}. 
#' @param X HiC experiment.
#' @param parallel Run in parallel? Must register parallel backend first. 
#' @export

diagNormalize = function(X,parallel=F){
  # Get scale factors
  S = scaleFactors(X=X,parallel=parallel);
  # Apply normalization
  Z = apply.DCSexp(X=X,f=diagScale,S=S,parallel=parallel);
  return(Z);
}

########################
# 2. Plot normalization
########################

#' Plot HiC Normalization Curves
#' 
#' Plots mean signal intensity by pixel separation before and after
#' HiC normalization. 
#' 
#' @param X Initial HiC experiment.
#' @param Z HiC Experiment after Normalization.
#' @param parallel Run in parallel? Must register parallel backend first. 
#' 
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @importFrom stats aggregate
#' @export

plotNormCurves = function(X,Z,parallel=F){
  # Before normalization
  Z.0 = apply.DCSexp(X=X,f=diagMeans,exp.out=F,parallel=parallel);
  Z.0 = do.call(rbind,Z.0);
  d = Sample = y = NULL;
  Z.0 = aggregate(y ~ Sample + d, data=Z.0, FUN=mean);
  Z.0$d = as.numeric(as.character(Z.0$d));
  
  q = ggplot(data=Z.0,aes(x=d,y=y,color=Sample)) + geom_line() +
    geom_point(size=1.5);
  q = q + theme_bw() + labs("x"="Distance between Loci",y="Mean Signal");
  q = q + scale_x_continuous(breaks=seq(from=0,to=max(Z.0$d),by=2))
  # y axis breaks
  b = 2^(seq(from=0,to=ceiling(log2(max(Z.0$y))),length.out=10));
  q = q + scale_y_continuous(breaks=b,labels=round(b,digits=2),trans="log2",limits=c(1,max(b)));
  q = q + theme(axis.title=element_text(size=12),title=element_text(size=14))
  q1 = q + ggtitle("Before Normalization");
  
  # After normalization
  Z.1 = apply.DCSexp(X=Z,f=diagMeans,exp.out=F,parallel=parallel);
  Z.1 = do.call(rbind,Z.1);
  Z.1 = aggregate(y ~ Sample + d, data=Z.1, FUN=mean);
  Z.1$d = as.numeric(as.character(Z.1$d));
  
  q = ggplot(data=Z.1,aes(x=d,y=y,color=Sample)) + geom_line() +
    geom_point(size=1.5);
  q = q + theme_bw() + labs("x"="Distance between Loci",y="Mean Signal");
  q = q + scale_x_continuous(breaks=seq(from=0,to=max(Z.0$d),by=2));
  # y axis breaks
  q = q + scale_y_log10(breaks=b,labels=round(b,digits=2),limits=c(1,max(b)));
  q = q + theme(axis.title=element_text(size=12),title=element_text(size=14))
  q2 = q + ggtitle("After Normalization");
  # Plot
  Q = cowplot::plot_grid(plotlist=list(q1,q2),align="h",nrow=1);
  return(Q);
}

