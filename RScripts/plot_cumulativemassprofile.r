
library(ggplot2)
library(grid)
library(emdbook)
source("function_plot-Mr-credreg.r")
source("functions_deason.r")

source("script_set-up.r")

# final.chain = readRDS(file = "../../../../PaperIII-chains/ErratumChains/modelparschaincombined-GCdata-correctSunmotion-Deason-newhier-210000-iterations-thinning20")

# final combined chain using Vasiliev data
final.chainVasiliev = readRDS(file = "../Results/modelpars-3-chainscombined-GCdata_noSgrGCsnoCrater_Vasiliev_2018-09-24_2018-09-24-Deason-rmin-15-2018-10-04_2018-10-04_2018-10-04_2018-10-05_2018-10-05")


# final combined chain using Gaia + HSTPROMO
final.chainGaiaHSTPROMO = readRDS(file = "../Results/modelpars-3-chainscombined-GCdata_noSgrGCsnoCrater_GaiaHSTPROMO_2018-10-05_2018-10-05-Deason-rmin-15-2018-10-05_2018-10-05_2018-10-05_2018-10-05_2018-10-06")

# # final chain with Gaia + HSTPROMO, rmin = 20
# final.chain20 = readRDS(file = "../Results/modelpars-3-chainscombined-GCdata_noSgrGCsnoCrater_GaiaHSTPROMO_2018-10-05_2018-10-05-Deason-rmin-20-2018-10-08_2018-10-08_2018-10-08_2018-10-08_2018-10-09")
# 
# # final chain with Vasiliev, rmin = 20
# final.chain20 = readRDS(file = "../Results/modelpars-3-chainscombined-GCdata_noSgrGCsnoCrater_Vasiliev_2018-09-24_2018-09-24-Deason-rmin-20-2018-10-09_2018-10-09_2018-10-09_2018-10-09_2018-10-09")


# axis label size
labsize = 1.25
# probability levels wanted
probnums = c(0.5, 0.75, 0.95)
credcols = rev(alpha(colour = "black", alpha = probnums))


source("script_previousstudies.r")

# colcompare = rainbow(1, start = 0.7, end=0.8)
mylwd = 1.5
mycex = 1.2
fillcolour = alpha("black", c(0.3, 0.2, 0.1))

# mars = c(5,5,2,20) # B ,L, T, R  <---- USE THESE IF YOU WANT THE LEGEND
mars = c(5,5,2,2) # <---- USE THESE IF YOU DON'T WANT THE LEGEND

#-------------- make cumulative mass profile plot
pdf(file = paste("../Figures/plot_cumulativemassprofile-", filename, "_", Sys.Date(), ".pdf", sep=""), useDingbats = FALSE, height=6, width=8)

# par( mar=mars, oma = outermars, xpd=TRUE)  #xpd = TRUE makes drawing outside margins possible
par(mfrow=c(1,2))

y.range=c(0,1.6)
plot.credreg(chain = final.chainVasiliev, dat = mydata, rmax = 160, printpdf = FALSE, pdf.name=NULL,
             regions = c(0.5, 0.75, 0.95), xyrange=TRUE, xrange=c(0,150), yrange=y.range,
             solids = fillcolour, cex.axis=mycex, cex=mycex, cex.lab=mycex, mymar = mars,
            plotminmax=TRUE, makelegend=FALSE)

points(comparethese[,1:2], pch=pchcompare, col=colcompare, cex=mycex)
arrows(x0 = comparethese[,1], x1 = comparethese[,1], y0 = comparethese[,2], y1=comparethese[,3], angle = 90,length = 0.1, col=colcompare, lwd=mylwd)
arrows(x0 = comparethese[,1], x1 = comparethese[,1], y0 = comparethese[,2], y1=comparethese[,4], angle = 90,length = 0.1, col=colcompare, lwd=mylwd)


par(mar = c(5,0,2,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("topleft", bg="white", legend = namecompare, pch=pchcompare, col=colcompare, title="Previous Studies", cex=0.8)



dev.off()



