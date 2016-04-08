##############################################################################
##############################################################################
# Example 3.13
# Particle Metropolis-Hastings for Swedish inflation/unemployment
# Reproduces Figure 3.9
#
# Copyright (c) 2016 Johan Dahlin [ johan.dahlin (at) liu.se ]
# Distributed under the MIT license.
#
##############################################################################
##############################################################################

#setwd("~/projects/phd-thesis/thesis-figures/ex-philips")

# Setup plot colors
library("RColorBrewer")
plotColors = brewer.pal(6, "Dark2");

# Settings for plotting
nMCMC   = 15000;
burnin  = 5000;


dpmh <- read.table("pmh0/philips_th_posterior.csv",header=T,sep=",");

cairo_pdf("ex-philips-pmh0.pdf",  height = 5, width = 8)
layout(matrix(1:4, 2, 2, byrow = TRUE))  
par(mar=c(4,5,1,1))

# Histograms with kernel density estimates
hist(dpmh$th0[burnin:nMCMC],breaks=floor(sqrt(nMCMC-burnin)),main="",freq=F,col=rgb(t(col2rgb(plotColors[1]))/256,alpha=0.25),border=NA,xlab=expression(phi),ylab="posterior estimate",xlim=c(0.5,1),ylim=c(0,8));
abline(v=mean(dpmh$th0[burnin:nMCMC]),lty="dotted");
lines(density(dpmh$th0[burnin:nMCMC],kernel="e",from=0.5,to=1),lwd=2,col=plotColors[1])

# Prior for phi
grid = seq(0.5,1,0.01);
dist = dnorm(grid,0.8,0.10)
lines(grid,dist,lwd=1,col="grey30")

hist(dpmh$th1[burnin:nMCMC],breaks=floor(sqrt(nMCMC-burnin)),main="",freq=F,col=rgb(t(col2rgb(plotColors[2]))/256,alpha=0.25),border=NA,xlab=expression(alpha),ylab="posterior estimate",xlim=c(0.0,1.0),ylim=c(0,3));
abline(v=mean(dpmh$th1[burnin:nMCMC]),lty="dotted");
lines(density(dpmh$th1[burnin:nMCMC],kernel="e",from=0.0,to=1.0),lwd=2,col=plotColors[2])

# Prior for alpha
grid = seq(0.0,1,0.01);
dist = dnorm(grid,0.5,0.2)
lines(grid,dist,lwd=1,col="grey30")

hist(dpmh$th2[burnin:nMCMC],breaks=floor(sqrt(nMCMC-burnin)),main="",freq=F,col=rgb(t(col2rgb(plotColors[3]))/256,alpha=0.25),border=NA,xlab=expression(beta),ylab="posterior estimate",xlim=c(-0.05,0.05),ylim=c(0,100));
abline(v=mean(dpmh$th2[burnin:nMCMC]),lty="dotted");
lines(density(dpmh$th2[burnin:nMCMC],kernel="e",from=-0.04,to=0.04),lwd=2,col=plotColors[3])

# Prior for beta
grid = seq(-0.04,0.1,0.04);
dist = dnorm(grid,0,0.1)
lines(grid,dist,lwd=1,col="grey30")

hist(dpmh$th3[burnin:nMCMC],breaks=floor(sqrt(nMCMC-burnin)),main="",freq=F,col=rgb(t(col2rgb(plotColors[4]))/256,alpha=0.25),border=NA,xlab=expression(sigma[v]),ylab="posterior estimate",xlim=c(0.2,0.35),ylim=c(0,50));
abline(v=mean(dpmh$th3[burnin:nMCMC]),lty="dotted");
lines(density(dpmh$th3[burnin:nMCMC],kernel="e",from=0.2,to=0.35),lwd=2,col=plotColors[4])

# Prior for sigma_v
grid = seq(0.2,0.35,0.01);
dist = dgamma(grid,shape=2.0,rate=4.0)
lines(grid,dist,lwd=1,col="grey30")
dev.off()