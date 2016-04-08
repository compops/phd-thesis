##############################################################################
##############################################################################
# Example 4.1
# Particle Metropolis-Hastings for Swedish inflation/unemployment
# Reproduces Figure 4.5
#
# Copyright (c) 2016 Johan Dahlin [ johan.dahlin (at) liu.se ]
# Distributed under the MIT license.
#
##############################################################################
##############################################################################


#setwd("~/projects/phd-thesis/thesis-figures/ex-philips")

# Setup plot colors
library("RColorBrewer")
plotColors = brewer.pal(8, "Dark2");

# Settings for plotting
nMCMC   = 10000;
burnin  = 2500;
plotM   = seq(burnin,nMCMC,1)
toPlot  = seq(8000,8300,1)

# Set working directory and read the results from file
dpmh1 <- as.numeric( read.table("qpmh2/philips_post_qpmh2.csv",header=F,sep=",")[,3]   )
dpmh2 <- c( rep(0,2500), as.numeric( read.table("qpmh2/philips_post_qpmh2zv.csv",header=F,sep=",")[,3] ) )
dpmh3 <- as.numeric( read.table("pmh0/philips_th_posterior.csv",header=T,sep=",")$th2  )

cairo_pdf("ex-philips-qpmh2.pdf",  height = 9, width = 8)

layout(matrix(1:9, 3, 3, byrow = FALSE))  
par(mar=c(4,5,1,1)) 

# Trace plots and ACFs
plot(toPlot,dpmh1[toPlot],type="l",ylab=expression(beta),xlab="iteration",col=plotColors[5],bty="n",ylim=c(-0.05,0.05));
polygon(c(toPlot,rev(toPlot)),c(dpmh1[toPlot],rep(-0.050,length(toPlot))),border=NA,col=rgb(t(col2rgb(plotColors[5]))/256,alpha=0.25))

plot(toPlot,dpmh2[toPlot],type="l",ylab=expression(beta),xlab="iteration",col=plotColors[7],bty="n",ylim=c(-0.05,0.05));
polygon(c(toPlot,rev(toPlot)),c(dpmh2[toPlot],rep(-0.050,length(toPlot))),border=NA,col=rgb(t(col2rgb(plotColors[7]))/256,alpha=0.25))

plot(toPlot,dpmh3[toPlot],type="l",ylab=expression(beta),xlab="iteration",col=plotColors[8],bty="n",ylim=c(-0.05,0.05))
polygon(c(toPlot,rev(toPlot)),c(dpmh3[toPlot],rep(-0.050,length(toPlot))),border=NA,col=rgb(t(col2rgb(plotColors[8]))/256,alpha=0.25))


foo = acf(dpmh1[burnin:nMCMC],lag.max=120,plot=FALSE);
plot(foo$lag,foo$acf,type="l",ylab=expression("acf of "* beta),col=plotColors[5],lwd=1,xlab="lag",bty="n",ylim=c(-0.2,1))
polygon(c(foo$lag,rev(foo$lag)),c(foo$acf,rep(-0.20,length(foo$lag))),border=NA,col=rgb(t(col2rgb(plotColors[5]))/256,alpha=0.25))
abline(h=-1.96/sqrt(length(dpmh1[burnin:nMCMC])),lty="dotted")
abline(h=1.96/sqrt(length(dpmh1[burnin:nMCMC])),lty="dotted")

foo = acf(dpmh2[burnin:nMCMC],lag.max=120,plot=FALSE);
plot(foo$lag,foo$acf,type="l",ylab=expression("acf of "* beta),col=plotColors[7],lwd=1,xlab="lag",bty="n",ylim=c(-0.2,1))
polygon(c(foo$lag,rev(foo$lag)),c(foo$acf,rep(-0.20,length(foo$lag))),border=NA,col=rgb(t(col2rgb(plotColors[7]))/256,alpha=0.25))
abline(h=-1.96/sqrt(length(dpmh1[burnin:nMCMC])),lty="dotted")
abline(h=1.96/sqrt(length(dpmh1[burnin:nMCMC])),lty="dotted")

foo = acf(dpmh3[burnin:nMCMC],lag.max=120,plot=FALSE);
plot(foo$lag,foo$acf,type="l",ylab=expression("acf of "* beta),col=plotColors[8],lwd=1,xlab="lag",bty="n",ylim=c(-0.2,1))
polygon(c(foo$lag,rev(foo$lag)),c(foo$acf,rep(-0.20,length(foo$lag))),border=NA,col=rgb(t(col2rgb(plotColors[8]))/256,alpha=0.25))
abline(h=-1.96/sqrt(length(dpmh1[burnin:nMCMC])),lty="dotted")
abline(h=1.96/sqrt(length(dpmh1[burnin:nMCMC])),lty="dotted")


# Histograms with kernel density estimates
grid = seq(-0.05,0.05,0.05);
dist = dnorm(grid,0,0.05)

hist(dpmh1[burnin:nMCMC],breaks=floor(sqrt(nMCMC-burnin)),main="",freq=F,col=rgb(t(col2rgb(plotColors[5]))/256,alpha=0.25),border=NA,xlab=expression(beta),ylab="posterior estimate",xlim=c(-0.05,0.05),ylim=c(0,90))
abline(v=mean(dpmh1[burnin:nMCMC]),lty="dotted");
lines(density(dpmh1[burnin:nMCMC],kernel="e",from=-0.05,to=0.05),lwd=2,col=plotColors[5])
lines(grid,dist,lwd=1,col="grey30")
text(0.050,90,pos=2,"qPMH2")

hist(dpmh2[burnin:nMCMC],breaks=floor(sqrt(nMCMC-burnin)),main="",freq=F,col=rgb(t(col2rgb(plotColors[7]))/256,alpha=0.25),border=NA,xlab=expression(beta),ylab="posterior estimate",xlim=c(-0.05,0.05),ylim=c(0,90))
abline(v=mean(dpmh2[burnin:nMCMC]),lty="dotted");
lines(density(dpmh2[burnin:nMCMC],kernel="e",from=-0.05,to=0.05),lwd=2,col=plotColors[7])
lines(grid,dist,lwd=1,col="grey30")
text(0.050,90,pos=2,"qPMH2-ZV")

hist(dpmh3[burnin:nMCMC],breaks=floor(sqrt(nMCMC-burnin)),main="",freq=F,col=rgb(t(col2rgb(plotColors[8]))/256,alpha=0.25),border=NA,xlab=expression(beta),ylab="posterior estimate",xlim=c(-0.05,0.05),ylim=c(0,90))
abline(v=mean(dpmh3[burnin:nMCMC]),lty="dotted");
lines(density(dpmh3[burnin:nMCMC],kernel="e",from=-0.05,to=0.05),lwd=2,col=plotColors[8])
lines(grid,dist,lwd=1,col="grey30")
text(0.050,90,pos=2,"PMH0")

dev.off()