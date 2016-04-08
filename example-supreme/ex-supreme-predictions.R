##############################################################################
##############################################################################
# Example 5.2
# US Supreme court ideological leaning
# Reproduces Figure 5.2
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

thetahat <- as.matrix(posterior1)[,1:9]
alphahat <- as.matrix(posterior1)[,seq(10,95,2)]
betahat <- as.matrix(posterior1)[,seq(11,95,2)]

nSim     <- 10000
alphaSim <- sample(alphahat,nSim,replace=TRUE)
betaSim  <- sample(betahat,nSim,replace=TRUE)
thetaSim <- apply( thetahat, 2, sample, nSim, replace=TRUE )

test = sweep( sweep(thetaSim, 1, betaSim, "*"), 1, -alphaSim, "+")

cairo_pdf("ex-supreme-predictions.pdf", height = 4, width = 8)

layout(matrix(1:2, 1, 2, byrow = TRUE))  
par(mar=c(4,5,1,1))

##
foo <- density( rowSums( test > 0 ) + 0.4 * rnorm(nSim), to=9)
hist( rowSums( test > 0 ), main="", col="darkgrey", border=NA, freq=FALSE, xlab="no. liberal votes", ylab="density", ylim=c(0,0.8), xlim=c(0,10) )
lines( foo, lwd=2, col=plotColors[1])

idx <- which( foo$x > 4.5)[1] -1 
polygon(c(foo$x[-(1:idx)],rev(foo$x[-(1:idx)])),c(foo$y[-(1:idx)],rep(0,length(foo$x[-(1:idx)]))),
        border=NA,col=rgb(t(col2rgb(plotColors[1]))/256,alpha=0.5))

# Probability of liberal majority (0.4248)
mean( rowSums( test > 0 ) > 4 )

## Scalia dies, What happends if someone more liberal replaces him?

thetaSim2     <- thetaSim
thetaSim2[,1] <- thetaSim[,1] - 1
test2 = sweep( sweep(thetaSim2, 1, betaSim, "*"), 1, -alphaSim, "+")

foo <- density( rowSums( test2 > 0 ) + 0.4 * rnorm(nSim), to=9)
hist( rowSums( test2 > 0 ), main="", col="darkgrey", border=NA, freq=FALSE, xlab="no. liberal votes", ylab="density", ylim=c(0,0.8), xlim=c(0,10) )
lines( foo, lwd=2, col=plotColors[2])

idx <- which( foo$x > 4.5)[1] -1 
polygon(c(foo$x[-(1:idx)],rev(foo$x[-(1:idx)])),c(foo$y[-(1:idx)],rep(0,length(foo$x[-(1:idx)]))),
        border=NA,col=rgb(t(col2rgb(plotColors[2]))/256,alpha=0.5))


# Probability of liberal majority (0.4654)
mean( rowSums( test2 > 0 ) > 4 )

dev.off()