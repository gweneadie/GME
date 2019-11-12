# plot prior distribution for alpha, given hyperprior parameters determined from other data

# source hyperprior set-up file
source("script_set-hyperpriors.r")

# function of x for plotting alpha prior with function curve()
plotalpharior = function(x, ppars){
  
  rmin = ppars[1]
  myc = ppars[2]
  p = ppars[3]
  n = ppars[4]
  nco = ppars[5]
  
  dgamma(x=(x-3), shape = (n+myc), scale=1/(nco+p))
  
}

#uncomment if you want to save plot
pdf("../Figures/priordist_alpha.pdf", width = 5, height=5, useDingbats = FALSE)

# set margins if you want
par(mar=c(5,5,2,1)) #bottom, left, top, right

curve( plotalpharior(x, ppars=apriorpars), ylab=expression(p(alpha~"|"~"b,c,p")), xlab=expression(alpha), lwd=2, cex.lab=2, cex.axis=2, ylim=c(0,10), xlim=c(3,5))

grid()

legend("topright", c(expression(b==0.4~"kpc"), expression(c==0.001), expression(p==0.001), expression(r[min]==~bquote(rmin))), lty=c(NA,NA,NA), lwd=c(NA,NA,NA), col=c(NA,NA,NA), cex=1.5)

# to turn off the device, or to close the pdf file
dev.off()

