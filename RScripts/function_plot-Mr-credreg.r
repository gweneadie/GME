# function to plot credible regions
# assumes parameters in chain are in regular space
library("emdbook")

plot.credreg <- function( chain, dat, rmax, pdf.name, modelM.r=M.r, rmin=1e-3, n.r=1e2, npars=2, printpdf=TRUE,
                          regions=c(0.5, 0.75, 0.95), plottrue=FALSE, plotminmax=FALSE, xyrange=FALSE,
                          xrange=FALSE, yrange=FALSE, r.true=FALSE, M.true=FALSE, othersim=FALSE, truepars=NULL,
                          # default colours for credible region bounds
                          solids = c("cadetblue4", "cadetblue3", "cadetblue1"),
                          ylab=expression( M(r<R)~(10^12~M[sol])), xlab="R (kpc)", 
                          mymar = c(5,7,0.5,5), makelegend=TRUE, ...){

  # r.values
  r.values = lseq(rmin, rmax, length.out=n.r)
  
  # credible region bounds
  drop.perc = (1 - regions)/2

  # columns = r , rows = M(r) in units of 10^12 solar masses
  Mr.values = 2.325e-3*sapply(X=r.values, FUN=modelM.r, pars= chain[ ,1:npars] )
  colnames(Mr.values) = paste( "r", seq(1,length(r.values)), sep="")
  
  # sort the M(r) values
  sorted.Mr = apply(X=Mr.values, MARGIN=2, FUN=sort)
  
  lower.creds = sapply( X=drop.perc, FUN=function(x) ( sorted.Mr[ -( 1:(x*nrow(sorted.Mr)) ) , ] )[1, ] )
  colnames(lower.creds) = as.character(100*regions)
  
  upper.creds = apply( X=cbind( drop.perc, regions ), MARGIN=1,
                       FUN=function(x) ( sorted.Mr[ -( 1:( sum(x) * nrow( sorted.Mr ) ) ), ] )[1, ] )
  colnames(upper.creds) = colnames(lower.creds)
  
  
  # r values for drawing sigma bounds (pass to polygon)
  poly.r = c( r.values, rev(r.values) )
  
  # credible bound M(r) values for plotting (pass to polygon)
  polyM = rbind( lower.creds, upper.creds[ nrow(upper.creds):1 , ] )
  
  # r values for plotting
  rs = c(r.values, rev(r.values))
  Ms = rep(0, n.r)
  
  # plot ranges
  if( xyrange ){ xrange = xrange ; yrange = yrange}else{
    xrange = c(rmin,rmax)
    yrange = c(0,max(Mr.values))
  }

if(printpdf){  pdf( pdf.name, useDingbats=FALSE  ) }

  par(mar=mymar)
  
  p = plot( r.values, (Ms), type="n", xlim=xrange, ylim=yrange,
            xlab=xlab, ylab=ylab, cex.lab=1.8, ...)
  #add a background grid
  grid()

  # add credible regions
  for( j in length(regions):1 ){
    polygon(x=rs, y=polyM[,j], col=solids[j], border=solids[j]) 
  }
  grid()
  
  if( makelegend ){
    # add legend without rmin/rmax
    legend( "topleft",  paste(colnames(lower.creds), "%" ), bg="white",
            lwd=c(8,8,8), lty=c(1,1,1), col=solids,  cex=0.8)
    
    # use if want to plot min/max of data position and key in legend
    if( plotminmax ){ 
      abline(v=max(dat[,"Rgc"]), lty=2)
      abline(v=min(dat[,"Rgc"]), lty=2)
      #covers previous legend
      legend( "topleft",  c(paste(colnames(lower.creds), "%" ), expression(r[min/max~data])), bg="white",
              lwd=c(8,8,8,1.3), lty=c(1,1,1,2), col=c(solids,"black"),  cex=0.8)
      
    }
  }


  # add true M(r) if needed
  if( plottrue & !othersim ){
    lines(x=r.values, y=2.325e-3*sapply(r.values, FUN=M.r, pars=truepars), col="red", cex=1.3) 
  }
  if( plottrue & othersim ){
    lines(x=r.true, y=M.true, col="red", cex=1.3)
  }
  

  print(p)
  
if(printpdf){  dev.off() }
  
}