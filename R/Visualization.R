########################
# Component Functions
########################

# Conventions:
# 1. Pixel coordinates are formatted as (i,j), where i is row and j is column of
# the contact matrix
# 2. The contact matrix is taken as upper-triangular, i.e. i<= j

#' Heatmap Pixel Coords
#' 
#' Given a sequence of foci, enumerates coordinates of all pixels included in
#' the heatmap.
#' @param s Sequence of foci.
#' @return Data frame containing pixel coordinates (i,j) and unique identifiers. 
#' @import plyr
Pixels = function(s){
  # Observed sequence
  n = length(s);
  l = min(s);
  u = max(s);
  # Coord i
  i = unlist(plyr::llply(.data=s,.fun=function(x){rep(x,times=(n-x+l))}));
  # Coord j
  j = unlist(plyr::llply(.data=s,.fun=function(x){seq(from=x,to=u)}));
  # Assigns labels to each pixel of the parallelogram
  TriIndex = function(n){
    # Sequence
    lab = seq(1:n);
    # Index i
    i = unlist(plyr::llply(.data=lab,.fun=function(x){rep(x,times=(n-x+1))}));
    # Index j
    j = unlist(plyr::llply(.data=lab,.fun=function(x){seq(from=x,to=n)}))
    # Combine
    R = apply(cbind(i,j),MARGIN=1,function(x){paste0(x[1],".",x[2])});
    return(R);
  }
  # Add index in parallelogram
  id = TriIndex(n);
  P = data.frame("i"=i,"j"=j,"id"=id);
  return(P);
}

#' X-Coord of Pixel Vertices
#' 
#' Starting from the left-most vertex and proceeding clockwise, enumerates the
#' x-coordinates of a pixel's vertices. 
#' @param r Pixel coords (i,j);
polyX = function(r){
  d = r[2]-r[1];
  x1 = (d/2)+(r[1]-1);
  x2 = x1 + 1/2;
  x3 = x2 + 1/2;
  x4 = x3 - 1/2;
  return(c(x1,x2,x3,x4));
}

#' Y-Coord of Pixel Vertices
#' 
#' Starting from the left-most vertex and proceeding clockwise, enumerates the
#' y-coordinates of a pixel's vertices. 
#' @param r Pixel coords (i,j);
polyY = function(r){
  y1 = (r[2]-r[1]);
  y2 = y1 + 1;
  y3 = y2 - 1;
  y4 = y3 - 1;
  return(c(y1,y2,y3,y4));
}

#' Enumerate polygon vertices
#' 
#' Given the base height and flank, enumerates the vertices of pixels that
#' are connected to form the heatmap. 
#' @param h.u Length of upstream domain in pixels
#' @param h.d Length of downstream domain in pixels
#' @param l.u Length of upstream flank in pixels
#' @param l.d Length of downstream flank in pixels
#' @return List containing 1. the polygon centers, and 2. the polygon vertices
polyCoord = function(h.u,h.d,l.u,l.d){
  # Length of upstream domain in pixels
  U = h.u + l.u;
  # Length of downstream domain in pixels
  D = h.d + l.d;
  # Enumerate pixels 
  P = Pixels(seq(from=1,to=U+1+D));
  id = P[,c("id")];
  R = P[,c("i","j")];
  # Obtain vertex coordinates
  X = apply(R,MARGIN=1,FUN=polyX);
  Y = apply(R,MARGIN=1,FUN=polyY);
  # Take means to obtain centers
  cX = apply(X,MARGIN=2,FUN=mean);
  cY = apply(Y,MARGIN=2,FUN=mean);
  # Grid
  Grid = data.frame("x"=c(X),"y"=c(Y),"id" = rep(id,each=4),"pOrder"=rep(seq(1:4),times=nrow(R)));
  # Centers
  Cents = data.frame("x"=cX,"y"=cY,"id"=id);
  # Filter positions below and above limits
  Grid = Grid[!(Grid$y<0),];
  Out = list("Centers"=Cents,"Vertices"=Grid);
  return(Out);
}

#' Subset Flank Pixels
#' 
#' Given 1. coordinates for all pixels in a region, and 2. coordinates for non-flank
#' pixels in the region, subsets the coordinates of flank pixels.
#' @param R All pixels in region
#' @param S Non-flank pixels in region
subFlank = function(R,S){
  # Create identifiers
  aux = function(r){paste0(r[1],".",r[2])};
  id.r = apply(X=R,MARGIN=1,FUN=aux);
  id.s = apply(X=S,MARGIN=1,FUN=aux);
  # Flank positions
  is.f = !(id.r %in% id.s);
  # Flank coordinates
  Flank = R[is.f,];
  row.names(Flank) = NULL;
  return(Flank);
}

#' ggplot HiC Heatmap
#' 
#' Contains the ggplot commands for creating the HiC heatmap.
#' @param df Plotting frame
#' @param Foci Sequence of foci plotted
#' @param h.u Length of upstream domain in pixels
#' @param h.d Length of downstream domain in pixels
#' @param l.u Length of upstream flank in pixels

ggHiC = function(df,Foci,h.u,h.d,l.u){
  # X axis breaks
  n = length(Foci);
  xb = round(seq(from=1,to=n,length.out=10));
  xb = sort(unique(c(xb,h.u+l.u+1)));
  # X break labels
  xlab = Foci[xb];
  # Adjust label position
  xb = xb - 0.5;
  # ggplot
  q = ggplot(data=df,aes(x=x,y=y)) + geom_polygon(aes(fill=Signal,group=id));
  q = q + scale_x_continuous(name="Position",breaks=xb,labels=xlab);
  q = q + scale_y_continuous(name="",breaks=NULL);
  q = q + theme_bw() + theme(strip.background=element_rect(fill="white"),
                             axis.text.x = element_text(angle = 90,vjust=0.5),
                             axis.title=element_text(size=12),title=element_text(size=14));
  # Demarcate inter-domain region
  q = q + geom_abline(intercept=2*((h.u+l.u)),slope=-2,linetype="dashed",color="black",alpha=0.8,size=1.2);
  q = q + geom_abline(intercept=-2*((h.u+l.u)+1),slope=2,linetype="dashed",color="black",alpha=0.8,size=1.2);
  # Demarcate domains
  q = q + geom_abline(intercept=-2*l.u,slope=2,linetype="dashed",color="black",alpha=0.8,size=1.2);
  q = q + geom_abline(intercept=2*(h.u+l.u+1+h.d),slope=-2,linetype="dashed",color="black",alpha=0.8,size=1.2);
  return(q); 
}

########################
# Base Plotting Function
########################

#' Basic HiC Plot
#'
#' Creates a HiC heatmap from data of the form (i,j,stat), where stat is a continuous
#' statistic. 
#' 
#' @param A Chromosome-level data as (i,j,stat)
#' @param f Focus.
#' @param h.u Length of upstream domain in pixels
#' @param h.d Length of downstream domain in pixels
#' @param l.u Lenght of upstream flank in pixels
#' @param l.d Length of downstream flank in pixels
#' @import doMC ggplot2
#' @export

plotHicBasic = function(A,f,h.u,h.d,l.u,l.d){
  ## 1. Foci for inclusion in heatmap
  # Upstream
  U = f-(h.u+l.u);
  # Downstream
  D = f+(h.d+l.d);
  # Foci
  Foci = seq(from=U,to=D);
  n = length(Foci);
  ## 2. Coords of pixels in heatmap
  P = Pixels(Foci);
  # Merge in data
  df = merge(x=P,y=A,by=c("i","j"));
  colnames(df)[4]=c("Signal");
  # Signal Limits
  minSignal = min(df$Signal);
  maxSignal = max(df$Signal);
  ## 3. Enumerate vertices of pixels in heatmap
  V = polyCoord(h.u,h.d,l.u,l.d);
  Vertices = V[["Vertices"]];
  ## 4. Plotting
  # Merge HiC data with vertices
  df.p = merge(x=Vertices,y=df,by="id");
  # Sort to ensure pixel vertices are connected in the correct order
  df.p = df.p[order(df.p$id,df.p$pOrder),]
  # Create plot
  q = ggHiC(df=df.p,Foci=Foci,h.u=h.u,h.d=h.d,l.u=l.u);
  # Output
  Q = list("Plot"=q,"minSignal"=minSignal,"maxSignal"=maxSignal);
  return(Q);
}

########################
# Smoothed HiC Plot
########################

#' Loess Smooth against Separation
#' 
#' @param X Matrix of the form (i,j,Signal)
#' @return Smoothed matrix.
Smooth = function(X){
  if(nrow(X)==1){
    return(X)
  } else {
    Out = X;
    # Calculate distance
    X$d = abs(X$i-X$j);
    # Regression
    M = loess(Signal~d,data=X);
    # Fitted values
    Out$Signal = fitted(M);
    return(Out);
  }
}

#' HiC Plot with Domain Smoothing
#' 
#' @param A Chromosome-level data as (i,j,stat)
#' @param f Focus.
#' @param h.u Length of upstream domain in pixels
#' @param h.d Length of downstream domain in pixels
#' @param l.u Lenght of upstream flank in pixels
#' @param l.d Length of downstream flank in pixels
#' @import doMC ggplot2
#' @export

plotHicSmoothed = function(A,f,h.u,h.d,l.u,l.d){
  # Relabel
  colnames(A) = c("i","j","Signal");
  ## 1. Plotting domains
  # Overall plotting domain
  a = f - (h.u+l.u);
  b = f + (h.d+l.d);
  Foci = seq(from=a,to=b);
  P = Pixels(Foci);
  R = P[,c(1:2)];
  # Domain subdivisions
  U = Pixels(seq(from=f-h.u,to=f-1))[,c(1:2)]; # Upstream
  D = Pixels(seq(from=f+1,to=f+h.d))[,c(1:2)]; # Downstream
  I = expand.grid(seq(from=f-h.u,to=f),seq(from=f,to=f+h.d));
  colnames(I) = c("i","j");
  # Non-flank pixels
  S = rbind(U,D,I);
  # Flank pixels
  Flank = subFlank(R=R,S=S);
  # Check if flank exists
  flankExists = (nrow(Flank)>0);
  ## 2. Merge in data
  U.data = merge(x=U,y=A,by=c("i","j"));
  D.data = merge(x=D,y=A,by=c("i","j"));
  I.data = merge(x=I,y=A,by=c("i","j"));
  if(flankExists){F.data = merge(x=Flank,y=A,by=c("i","j"));}
  ## 3. Smooth data
  # Spline smoothing
  U.data = Smooth(X=U.data);
  D.data = Smooth(X=D.data);
  I.data = Smooth(X=I.data);
  if(flankExists){F.data = Smooth(X=F.data);}
  ## 4. Merge in data
  P1 = merge(x=P,y=U.data,by=c("i","j"));
  P2 = merge(x=P,y=D.data,by=c("i","j"));
  P3 = merge(x=P,y=I.data,by=c("i","j"));
  df = rbind(P1,P2,P3);
  if(flankExists){
    P4 = merge(x=P,y=F.data,by=c("i","j"));
    df = rbind(df,P4);
  }
  # Signal Limits
  minSignal = min(df$Signal);
  maxSignal = max(df$Signal);
  # 3. Enumerate vertices of pixels in heatmap
  V = polyCoord(h.u,h.d,l.u,l.d);
  Vertices = V[["Vertices"]];
  ## 4. Plotting
  # Merge HiC data into lattice
  df.p = merge(x=Vertices,y=df,by="id");
  # Sort to ensure pixel vertices are connected in the correct order
  df.p = df.p[order(df.p$id,df.p$pOrder),]
  # Create plot
  q = ggHiC(df=df.p,Foci=Foci,h.u=h.u,h.d=h.d,l.u=l.u);
  # Output
  Q = list("Plot"=q,"minSignal"=minSignal,"maxSignal"=maxSignal);
  return(Q);
}

########################
# Grouped HiC Plot
########################

#' HiC Plot by Group
#' 
#' Creates a HiC heatmap faceted by group.
#' @param A Chromosome-level data as  (i,j, group1 stat, group2 stat, etc.)
#' @param f Focus.
#' @param h.u Length of upstream domain in pixels
#' @param h.d Length of downstream domain in pixels
#' @param l.u Lenght of upstream flank in pixels
#' @param l.d Length of downstream flank in pixels
#' @param Smooth Should pixel intensity be smoothed?
#' @import doMC ggplot2
#' @export

plotHicGrouped = function(A,f,h.u,h.d,l.u,l.d,smooth=F){
  # Groups
  L = sort(names(A[,3:ncol(A)]));
  ng = length(L);
  # Loop over groups
  G = foreach(i=1:ng) %do% {
    # Subset group columns
    Sub = A[,c("i","j",L[i])];
    if(smooth){
      Q = plotHicSmoothed(A=Sub,f=f,h.u=h.u,h.d=h.d,l.u=l.u,l.d=l.d);
    } else {
      Q = plotHicBasic(A=Sub,f=f,h.u=h.u,h.d=h.d,l.u=l.u,l.d=l.d);
    }
    return(Q);
  } # End loop over groups
  
  ## Harmonize scales
  signalScales = foreach(i=1:ng,.combine=rbind) %do% {
    return(c(G[[i]]$minSignal,G[[i]]$maxSignal))
  }
  # Min Signal
  minSignal = min(signalScales[,1]);
  minSignal = floor(minSignal*100)/100;
  # Max Signal
  maxSignal = max(signalScales[,2]);
  maxSignal = ceiling(maxSignal*100)/100;
  
  plotList = foreach(i=1:ng) %do% {
    q = G[[i]]$Plot + scale_fill_distiller(palette="Spectral",limits=c(minSignal,maxSignal));
    q = q + ggtitle(L[i]);
    return(q)
  }
  
  # Arrange into grid
  qout = cowplot::plot_grid(plotlist=plotList,ncol=1);
  return(qout);
}
