##############################################################################
##############################################################################
# Example 3.6
# Particle filtering for Swedish inflation/unemployment
# Reproduces Figure 3.7
#
# Copyright (c) 2016 Johan Dahlin [ johan.dahlin (at) liu.se ]
# Distributed under the MIT license.
#
##############################################################################
##############################################################################

setwd("~/projects/phd-thesis/thesis-figures/ex-philips")

# Setup plot colors
library("RColorBrewer")
plotColors = brewer.pal(6, "Dark2");

xhatm <- read.table("smc/philips_smc_llgrid.csv",header=F,sep=",",stringsAsFactors=F)[-1,]
grid  <- seq(0.01,0.99,0.01)

cairo_pdf("ex-philips-llscore.pdf",  height = 8, width = 8)

layout(matrix(1:2, 2, 1, byrow = TRUE))  
par(mar=c(4,5,1,1))

plot(grid[1:92],xhatm[1:92,1],type="l",xlab=expression(phi),ylab="log-likelihood",col=plotColors[3],bty="n",ylim=c(-50,-45),xlim=c(0,1))
polygon(c(grid[1:92],rev(grid[1:92])),c(xhatm[1:92,1],rep(-50,length(grid[1:92]))),border=NA,col=rgb(t(col2rgb(plotColors[3]))/256,alpha=0.25))
lines(grid[1:92],xhatm[1:92,1],lwd=1,col=plotColors[3])
abline(v=0.75,lty="dotted")

plot(grid,xhatm[,2],type="l",xlab=expression(phi),ylab=expression("gradient of log-posterior wrt " *phi) ,col=plotColors[4],bty="n",ylim=c(-400,400))
polygon(c(grid,rev(grid)),c(xhatm[,2],rep(-400,length(grid))),border=NA,col=rgb(t(col2rgb(plotColors[4]))/256,alpha=0.25))
lines(grid,xhatm[,2],lwd=1,col=plotColors[4])
abline(v=0.81,lty="dotted")
abline(h=0,lty="dotted")

dev.off()