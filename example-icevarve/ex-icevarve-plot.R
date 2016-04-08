##############################################################################
##############################################################################
# Example 1.1
# Ice varve thickness
# Reproduces Figure 1.1
#
# Copyright (c) 2016 Johan Dahlin [ johan.dahlin (at) liu.se ]
# Distributed under the MIT license.
#
##############################################################################
##############################################################################
# Setup plot colors
library("RColorBrewer")
plotColors = brewer.pal(6, "Dark2");

#setwd("~/projects/phd-thesis/thesis-figures/ex-icevarve")

d          = read.table("icevarve.txt",sep=",",header=F,stringsAsFactors=F)
postGP     = read.table("ex-icevarve-posterior.csv",sep=",",header=T)[,-1]
postPMH    = read.table("icevarve_pmh_xhats.csv",sep=",",header=F)$V1
postPMHvar = read.table("icevarve_pmh_xhats_var.csv",sep=",",header=F)$V1

grid  = rev(seq(9883,9250))

cairo_pdf("ex-icevarve.pdf",  height = 10, width = 8)

layout(matrix(1:3, 3, 1, byrow = TRUE))  
par(mar=c(4,5,1,1))

plot(grid,d$V1,type="l",col=plotColors[1],bty="n",ylim=c(0,150),xlab="year (BC)",ylab="ice varve thickness",xlim=c(9900,9200),xaxt="n")
polygon( c(grid,rev(grid)), c(d$V1, rep(0,length(grid))),border=NA,col=rgb(t(col2rgb(plotColors[1]))/256,alpha=0.25))
axis(1,at=rev(seq(9900,9200,-100)))

plot(grid,postGP$X0,type="l",bty="n",ylim=c(0,150),lwd=2,col=plotColors[2],xlab="year (BC)",ylab="predictive posterior",xlim=c(9900,9200),xaxt="n")
axis(1,at=rev(seq(9900,9200,-100)))
x2max = apply( cbind(rep(0,length(grid)),rev(postGP$X1)), 1, max )

polygon( c(grid,rev(grid)), c(postGP$X2, x2max),border=NA,col=rgb(t(col2rgb(plotColors[2]))/256,alpha=0.25))
points(grid,d$V1,pch=19,cex=0.25)
lines(grid,postGP$X0,lwd=2,col=plotColors[2])

grid2 <- grid[-(1:2)]
x2max = apply( cbind(rep(0,length(grid2)),rev(postPMH[-(1:2)]-1.96*sqrt(postPMHvar[-(1:2)]))), 1, max )
plot(grid2,postPMH[-(1:2)],type="l",bty="n",ylim=c(0,150),lwd=2,col=plotColors[3],xlab="year (BC)",ylab="predicted thickness",xlim=c(9900,9200),xaxt="n")
polygon( c(grid2,rev(grid2)), c(postPMH[-(1:2)]+1.96*sqrt(postPMHvar[-(1:2)]), x2max),border=NA,col=rgb(t(col2rgb(plotColors[3]))/256,alpha=0.25))
axis(1,at=rev(seq(9900,9200,-100)))
points(grid,d$V1,pch=19,cex=0.25)
lines(grid2,postPMH[-(1:2)],lwd=2,col=plotColors[3])

dev.off()