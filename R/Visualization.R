########################
# Component Functions
########################

#' Matrix Exclusion
#' 
#' Find all rows in matrix A NOT in matrix B
#' @param A First matrix.
#' @param B Second matrix.
#' @importFrom data.table data.table setkey

mExclusion = function(A,B){
  A = data.table::data.table(A);
  data.table::setkey(x=A);
  B = data.table::data.table(B);
  data.table::setkey(x=B);
  D = A[!B];
  return(D);
}

#' Heatmap Pixel Coords
#' 
#' Given a sequence of foci, enumerates coordinates of all pixels included in 
#' the heatmap.
#' @param s Sequence of foci.
#' @return Data.frame containing pixel coordinates (i,j) and a unique identifier in
#'   the form (1+i-min(i),1+j-min(i)).
#'   
#' @importFrom plyr llply
#' @importFrom data.table data.table setkey

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
  P = data.table::data.table("i"=i,"j"=j,"id"=id);
  data.table::setkey(P,i,j);
  return(P);
}

#' X-Coord of Pixel Vertices
#' 
#' Starting from the left-most vertex and proceeding counter-clockwise,
#' enumerates the x-coordinates of a pixel's vertices.
#' @param r Pixel coords (i,j);
#' @return Numeric vector of x-coordinates for heatmap.

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
#' Starting from the left-most vertex and proceeding counter-clockwise,
#' enumerates the y-coordinates of a pixel's vertices.
#' @param r Pixel coords (i,j);
#' @return Numeric vector of y-coordinates for heatmap.

polyY = function(r){
  y1 = (r[2]-r[1]);
  y2 = y1 + 1;
  y3 = y2 - 1;
  y4 = y3 - 1;
  return(c(y1,y2,y3,y4));
}

#' Enumerate polygon vertices
#' 
#' Given the base height and flank, enumerates the vertices of pixels that are
#' connected to form the heatmap.
#' 
#' @param h.u Length of upstream domain in pixels
#' @param h.d Length of downstream domain in pixels
#' @param l.u Length of upstream flank in pixels
#' @param l.d Length of downstream flank in pixels
#' @return List containing 1. the polygon centers, and 2. the polygon vertices. 
#'   \code{Centers} is a data.frame containing (x,y) coordinates of the pixel
#'   centers and the pixel id. \code{Vertices} is a data.frame containing (x,y)
#'   coordinates of the pixel vertices, the pixel id, and the order in which
#'   vertices are connected.
#'   
#' @importFrom data.table data.table setkey

polyCoord = function(h.u,h.d,l.u,l.d){
  # Length of upstream domain in pixels
  U = h.u + l.u;
  # Length of downstream domain in pixels
  D = h.d + l.d;
  # Enumerate pixels 
  P = Pixels(seq(from=1,to=U+1+D));
  id = c(P[,c("id")])[[1]];
  id = factor(id,levels=id,ordered=T);
  R = P[,c("i","j")];
  # Obtain vertex coordinates
  X = apply(R,MARGIN=1,FUN=polyX);
  Y = apply(R,MARGIN=1,FUN=polyY);
  # Take means to obtain centers
  cX = apply(X,MARGIN=2,FUN=mean);
  cY = apply(Y,MARGIN=2,FUN=mean);
  # Grid
  Grid = data.table::data.table("x"=c(X),"y"=c(Y),"id" = rep(id,each=4),"pOrder"=rep(seq(1:4),times=nrow(R)));
  data.table::setkey(Grid,id);
  # Centers
  Cents = data.table::data.table("x"=cX,"y"=cY,id);
  data.table::setkey(Cents,id);
  # Remove coordinates outside plot limits
  Grid = Grid[!(Grid$y<0),];
  Out = list("Centers"=Cents,"Vertices"=Grid);
  return(Out);
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
  x = y = Signal = id = NULL;
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
#' Creates a HiC heatmap displaying the value of the test statistic in each pixel.
#' 
#' @param B Data.table in the format (i,j,Signal).
#' @param f Plot focus, numeric.
#' @param h.u Length of upstream domain in pixels
#' @param h.d Length of downstream domain in pixels
#' @param l.u Lenght of upstream flank in pixels
#' @param l.d Length of downstream flank in pixels
#' @return List containing 1. the ggplot, 2. the minimal signal observed in the 
#'   plotting region, and 3. the maximal signal in the plotting region.
#'   
#' @import ggplot2

basicHicPlot = function(B,f,h.u,h.d,l.u,l.d){
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
  df = merge(x=P,y=B,by=c("i","j"));
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

#' Smooth Signal
#' 
#' Smooths the signal against pixel separation. If there are
#' fewer than 10 observations, a simple mean is taken. If between 10 and 20
#' observations, linear loess is applied. If greater than 20 observations,
#' quadratic loess is applied.
#' @param X Data table of the form (i,j,Signal)
#' @return Data table including the coordinates, original signal, pixel
#'   separation, and fitted values
#'   
#' @importFrom stats fitted loess

Smooth = function(X){
  # Calculate distance
  X$d = abs(X$i-X$j);
  if(nrow(X)<10){
    X$fit = mean(X$Signal);
  } else if (nrow(X)<20){
    # Regression
    M = loess(Signal~d,degree=1,data=X);
    # Fitted values
    X$fit = fitted(M);
  } else {
    # Regression
    M = loess(Signal~d,degree=2,data=X);
    # Fitted values
    X$fit = fitted(M);
  }
  return(X);
}

#' HiC Plot with Domain Smoothing
#' 
#' Creates a HiC heatmap displaying the fitted value of the test statistic smoothed
#' within domains by pixel separation. 
#' 
#' @param B Data.table in the format (i,j,Signal)
#' @param f Plot focus, numeric.
#' @param h.u Length of upstream domain in pixels
#' @param h.d Length of downstream domain in pixels
#' @param l.u Lenght of upstream flank in pixels
#' @param l.d Length of downstream flank in pixels
#' @return List containing 1. the ggplot, 2. the minimal signal observed in the 
#'   plotting region, and 3. the maximal signal in the plotting region.
#' 
#' @import ggplot2
#' @importFrom plyr llply
#' @importFrom data.table data.table setkey

smoothHicPlot = function(B,f,h.u,h.d,l.u,l.d){
  
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
  ## Domain subdivisions
  domains = list();
  # 1. Upper TAD
  if(h.u>0){D1 = Pixels(seq(from=f-h.u,to=f-1))[,c(1:2)]
    domains$D1 = D1;
  };
  # 2. Lower TAD
  if(h.d>0){D2 = Pixels(seq(from=f+1,to=f+h.d))[,c(1:2)]
    domains$D2 = D2;
  };
  # 3. Inter-TAD
  D3 = data.table::data.table(expand.grid(seq(from=f-h.u,to=f),seq(from=f,to=f+h.d)));
  colnames(D3) = c("i","j");
  data.table::setkey(D3);
  domains$D3 = D3;
  # 4. Upper Flank - Upper TAD
  if((l.u>0)&(h.u>0)){D4 = mExclusion(A=Pixels(seq(from=f-h.u-l.u,to=f-1))[,c(1:2)],B=D1)
    domains$D4 = D4;
  };
  # 5. Upper Flank - Lower TAD
  if((l.u>0)&(h.d>0)){D5 = data.table::data.table(expand.grid(seq(from=f-h.u-l.u,to=f-h.u-1),seq(from=f,to=f+h.d)))
    colnames(D5) = c("i","j");
    data.table::setkey(D5);
    domains$D5 = D5;
  };
  # 6. Upper Flank - Lower Flank
  if((l.u>0)&(l.d>0)){D6 = data.table::data.table(expand.grid(seq(from=f-h.u-l.u,to=f-h.u-1),seq(from=f+h.d+1,f+h.d+l.d)))
    colnames(D6) = c("i","j");
    data.table::setkey(D6);
    domains$D6 = D6;
  }
  # 7. Upper TAD - Lower Flank
  if((h.u>0)&(l.d>0)){D7 = data.table::data.table(expand.grid(seq(from=f-h.u,to=f),seq(from=f+h.d+1,f+h.d+l.d)))
    colnames(D7) = c("i","j");
    data.table::setkey(D7);
    domains$D7 = D7;
  };
  # 8. Lower TAD - Lower Flank
  if((h.d>0)&(l.d>0)){D8 = mExclusion(A=Pixels(seq(from=f+1,to=f+h.d+l.d))[,c(1:2)],B=D2)
    domains$D8 = D8;
  }
  
  ## 3. Merge in Data
  aux = function(x){
    return(merge(x=x,y=B));
  }
  domains = plyr::llply(.data=domains,.fun=aux);
  
  ## 4. Smooth Data
  domains = plyr::llply(.data=domains,.fun=Smooth);
  
  ## 5. Merge into overall data.table
  id = NULL;
  aux = function(x){
    merge(x=x,y=P);
  }
  domains = plyr::llply(.data=domains,.fun=aux);
  df = do.call(rbind,domains);
  data.table::setkey(df,id);
  
  # Signal Limits
  minSignal = min(df$fit);
  maxSignal = max(df$fit);
  
  ## 6. Enumerate vertices of pixels in heatmap
  V = polyCoord(h.u,h.d,l.u,l.d);
  Vertices = V[["Vertices"]];
  
  ## 4. Plotting
  # Merge HiC data into lattice
  df.p = merge(x=Vertices,y=df);
  # Sort to ensure pixel vertices are connected in the correct order
  df.p = df.p[order(df.p$id,df.p$pOrder),c("id","pOrder","x","y","fit")];
  colnames(df.p)[5] = "Signal";
  # Create plot
  q = ggHiC(df=df.p,Foci=Foci,h.u=h.u,h.d=h.d,l.u=l.u);
  # Output
  Q = list("Plot"=q,"minSignal"=minSignal,"maxSignal"=maxSignal);
  return(Q);
}

########################
# Plot HiC Statistic
########################

#' Plot HiC Statistic
#' 
#' Creates a heatmap displaying the value of a HiC statistic within each pixel. 
#' 
#' @param In A \code{DCSexp} or \code{DCSchr} object.
#' @param stat Name of statistic to plot. Should correspond with a column of
#'   \code{DCSchr@Stats}.
#' @param f Plot focus, numeric.
#' @param chr Chromosome, if providing a \code{DCSexp}.
#' @param h.u Length of upstream domain in pixels
#' @param h.d Length of downstream domain in pixels
#' @param l.u Lenght of upstream flank in pixels
#' @param l.d Length of downstream flank in pixels
#' @param smooth Smooth signal against pixel separation \eqn{d}?
#' @return List containing 1. the ggplot, 2. the minimal signal observed in the 
#'   plotting region, and 3. the maximal signal in the plotting region.
#' 
#' @import ggplot2
#' @importFrom plyr llply
#' @importFrom data.table data.table setkey
#' @export

plotHicStat = function(In,stat,f,chr,h.u,h.d,l.u=0,l.d=0,smooth=F){
  ## 0. Check input
  in.class = class(In);
  if(in.class=="DCSexp"){
    if(missing(chr)){stop("Must provide chromosome if input is a DCSexp.")};
    chr = as.character(chr);
    if(!(chr %in% In@chrs)){stop("Chromosome not present in the experiment.")};
    A = In@Data[[chr]]
  } else if(in.class=="DCSchr"){
    A = In;
  } else {stop("DCSexp or DCSchr expected as input.")};
  # Check focus
  if(!(f %in% A@Coord[,"i"])){stop("Focus not present on the chromosome.")};
  # Check statistic
  if(!(stat %in% colnames(A@Stats))){stop("Statistic not present for the chromosome.")};
  
  # Plotting frame
  i = j = NULL;
  B = data.table::data.table(A@Coord[,c("i","j")],"Signal"=A@Stats[,stat]);
  data.table::setkey(x=B,i,j);
  
  if(smooth){
    Q = smoothHicPlot(B=B,f=f,h.u=h.u,h.d=h.d,l.u=l.u,l.d=l.d)
  } else {
    Q = basicHicPlot(B=B,f=f,h.u=h.u,h.d=h.d,l.u=l.u,l.d=l.d);
  }
  return(Q);
}

########################
# Grouped HiC Plot
########################

#' Plot Log HiC Signal
#' 
#' Creates a heatmap displaying the log of the HiC signal in each pixel. The 
#' plots are facetted by group. For groups with multiple samples, the mean
#' signal is taken prior to the logarithm. 
#' 
#' @param In A \code{DCSexp} or \code{DCSchr} object.
#' @param f Plot focus, numeric.
#' @param chr Chromosome, if providing a \code{DCSexp}.
#' @param h.u Length of upstream domain in pixels
#' @param h.d Length of downstream domain in pixels
#' @param l.u Lenght of upstream flank in pixels
#' @param l.d Length of downstream flank in pixels
#' @param smooth Smooth signal against pixel separation \eqn{d}?
#' @return List containing 1. the group facetted heatmaps, 2. a list of the
#'   component plots, and 3. a data.frame containing the minimum and maximal 
#'   signal within the plotting region for each group.
#' 
#' @import ggplot2
#' @importFrom plyr llply
#' @importFrom data.table data.table setkey
#' @importFrom foreach "%do%" foreach 
#' @export

plotHicSignal = function(In,f,chr,h.u,h.d,l.u,l.d,smooth=F){
  ## 0. Check input
  in.class = class(In);
  if(in.class=="DCSexp"){
    if(missing(chr)){stop("Must provide chromosome if input is a DCSexp.")};
    chr = as.character(chr);
    if(!(chr %in% In@chrs)){stop("Chromosome not present in the experiment.")};
    A = In@Data[[chr]]
  } else if(in.class=="DCSchr"){
    A = In;
  } else {stop("DCSexp or DCSchr expected as input.")};
  # Check focus
  if(!(f %in% A@Coord[,"i"])){stop("Focus not present on the chromosome.")};
  
  ## 1. Groups
  L = sort(unique(A@groups));
  ng = length(L);
  # Function to summarize signal
  aux = function(x){log(mean(x)+1)};
  # Loop over groups
  i = j = NULL;
  G = foreach(i=1:ng) %do% {
    # Active group
    ag = L[i];
    # Subset samples in group, and take log of mean signal
    B = data.table::data.table(A@Coord[,c("i","j")],
                               "Signal"=plyr::aaply(.data=A@Counts[,A@groups %in% ag,drop=F],.margins=1,.fun=aux));
    data.table::setkey(x=B,i,j);
    if(smooth){
      Q = smoothHicPlot(B=B,f=f,h.u=h.u,h.d=h.d,l.u=l.u,l.d=l.d);
    } else {
      Q = basicHicPlot(B=B,f=f,h.u=h.u,h.d=h.d,l.u=l.u,l.d=l.d);
    }
    return(Q);
  } # End loop over groups
  
  ## Harmonize color scales
  signalScales = foreach(i=1:ng,.combine=rbind) %do% {
    return(c(G[[i]]$minSignal,G[[i]]$maxSignal))
  }
  signalScales = data.frame(signalScales);
  colnames(signalScales) = c("minSignal","maxSignal");
  rownames(signalScales) = NULL;
  signalScales$Group = L;
  
  # Min Signal
  minSignal = min(signalScales[,1]);
  minSignal = floor(minSignal*100)/100;
  # Max Signal
  maxSignal = max(signalScales[,2]);
  maxSignal = ceiling(maxSignal*100)/100;
  # Apply common color scale
  plotList = foreach(i=1:ng) %do% {
    q = G[[i]]$Plot + scale_fill_distiller(palette="Spectral",limits=c(minSignal,maxSignal));
    q = q + ggtitle(paste0("Group: ",L[i]));
    return(q)
  }
  # Arrange into grid
  qout = cowplot::plot_grid(plotlist=plotList,ncol=1);
  # Output
  Out = list("Plot"=qout,"plotList"=plotList,"signalScales"=signalScales);
  return(Out);
}
