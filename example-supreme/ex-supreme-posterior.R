##############################################################################
##############################################################################
# Example 3.8
# US Supreme court ideological leaning
# Reproduces Figure 3.8
#
# Copyright (c) 2016 Johan Dahlin [ johan.dahlin (at) liu.se ]
# Distributed under the MIT license.
#
##############################################################################
##############################################################################

#setwd("~/projects/phd-thesis/thesis-figures/ex-supreme")

# Setup plot colors
library("RColorBrewer")
plotColors = brewer.pal(8, "Dark2");
plotColors1 = c(plotColors,"#663333")
plotColors2 = c("#FFFFFF",plotColors,"#663333")

library("MCMCpack")

d <- as.matrix( read.table("supreme-court-2010-2015.csv",sep=",",stringsAsFactors=FALSE,header=TRUE) )
nCases <- dim(d)[1]

## US Supreme Court Example with inequality constraints

posterior1 <- MCMCirt1d(t(d),
                        theta.constraints=list(Scalia="+", Ginsburg="-"),
                        burnin=10000, mcmc=40000, thin=1, verbose=5000,
                        store.item=TRUE)
#summary(posterior1)

res <- as.matrix(posterior1)[,1:9]

cairo_pdf("ex-supreme-posteriors.pdf", height = 5, width = 8)
layout(matrix(1, 1, 1, byrow = TRUE))  
par(mar=c(4,5,1,1))

hist(res[,1],breaks=floor(sqrt(nMCMC-burnin)),main="",freq=F,
     col=rgb(t(col2rgb(plotColors1[1]))/256,alpha=0.25),border=NA,
     xlab="liberal/conservative score",ylab="density",xlim=c(-4,4),ylim=c(0,4))

lines(density(res[,1],kernel="e"),lwd=2,col=plotColors1[1])

for (ii in 2:9) {
  hist(res[,ii],breaks=floor(sqrt(nMCMC-burnin)),main="",freq=F,
       col=rgb(t(col2rgb(plotColors1[ii]))/256,alpha=0.25),border=NA,add=TRUE)
  lines(density(res[,ii],kernel="e"),lwd=2,col=plotColors1[ii])
}
legend(2,4,colnames(d),col=plotColors1,box.lwd=0,pch=15)
dev.off()

posteriorMeans <- colMeans(res) 
names(posteriorMeans) <- colnames(d)
sort(posteriorMeans)

upperLimit  = c(4,-1,4,0)
lowerLimit  = c(0,-4,-0.2,-2)

cairo_pdf("ex-supreme-traces.pdf", height = 5, width = 8)
layout(matrix(1:4, 2, 2, byrow = TRUE))  
par(mar=c(4,5,1,1))

for (ii in 3:4) {
  plot(res[1:1000,ii],type="l",
       col=plotColors1[ii-2],
       ylab=bquote(theta[.(ii-2)]), 
       xlab="iteration",bty='n',
       ylim=c(lowerLimit[ii-2],upperLimit[ii-2]))
  polygon(c(1:1000,rev(1:1000)),c(res[1:1000,ii],rep(lowerLimit[ii-2],1000)),
          border=NA,col=rgb(t(col2rgb(plotColors[ii-2]))/256,alpha=0.25))
  
  acfOut <- acf(res[1:1000,ii],plot=FALSE,lag.max=200)
  plot(acfOut$lag,acfOut$acf,type="l",
       col=plotColors1[ii-2],ylim=c(-0.2,1),
       xlab="lag",ylab=bquote("ACF of" ~ theta[.(ii-2)]),bty='n',lwd=2)
  polygon(c(acfOut$lag,rev(acfOut$lag)),c(acfOut$acf,rep(-0.2,length(acfOut$lag))),
          border=NA,col=rgb(t(col2rgb(plotColors[ii-2]))/256,alpha=0.25))
}
dev.off()