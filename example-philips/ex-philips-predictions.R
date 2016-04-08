##############################################################################
##############################################################################
# Example 5.1
# Particle Metropolis-Hastings for Swedish inflation/unemployment
# Reproduces Figure 5.1
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

d            <- read.table("philipscurve_sweden_1987_2015.csv",sep=",",header=T)
dates        <- as.Date( seq(as.POSIXct("1987-01-01 01:00:00 CET"), as.POSIXct("2015-12-01 01:00:00 CET"), by = "1 months") )
inflation    <- as.numeric( d$inflation    )
unemployment <- as.numeric( d$unemployment )

dpmh         <- read.table("pmh0/philips_th_posterior.csv",header=T,sep=",");
xhats        <- as.numeric( read.table("pmh0/philips_nairu_mean.csv",header=F,sep=",")$V1 );
xhatm        <- as.numeric( read.table("smc/state_bPF_N100.csv",header=T,sep=",")$xhatf );

posterior    <- dpmh[burnin:nMCMC,2:5]

future_unemployment <- unemployment[length(unemployment)] - 0.003 * seq(1,24,1)^2

nRuns <- 1000

un = matrix(xhatm[length(xhatm)], nrow=nRuns, ncol=24)
y  = matrix(inflation[length(inflation)], nrow=nRuns, ncol=24)

phi   = sample(posterior[,1], nRuns, replace=TRUE);
alpha = sample(posterior[,2], nRuns, replace=TRUE);
beta  = sample(posterior[,3], nRuns, replace=TRUE);
sigma = sample(posterior[,4], nRuns, replace=TRUE);

for ( ii in 1:nRuns ){
  for ( tt in 2:24 ) { 
    un[ii,tt] <- phi[ii] * un[ii,tt-1] + alpha[ii] / ( 1 + exp(-future_unemployment[tt-1])) + 1.0/(1.0+exp(abs(un[ii,tt-1]-future_unemployment[tt-1]))) * rnorm(1)
    y[ii,tt]  <- y[ii,tt-1] + beta[ii] * ( future_unemployment[tt] - un[ii,tt] ) + sigma[ii] * rnorm(1)
  }
}

dates        = seq(as.POSIXct("2015-01-01 01:00:00 CET"), as.POSIXct("2017-12-01 01:00:00 CET"), by = "1 months")
inflation    <- c( as.numeric( d$inflation    )[336:348], apply(y,2,mean)[-1] )
unemployment <- c( as.numeric( d$unemployment )[336:348], future_unemployment[-1] )

cairo_pdf("ex-philips-predictions.pdf",  height = 4, width = 8)

layout(matrix(c(1,2), 1, 2, byrow = TRUE))  
par(mar=c(4,5,1,1))

plot(as.Date(dates),unemployment,type="l",ylab="unemployment rate", xlab="date",col=plotColors[2],bty="n",xaxt="n", ylim=c(4,10), xlim=c(as.Date("2015-01-01 01:00:00 CET"),as.Date("2018-01-01 01:00:00 CET")) )

r         = as.POSIXct(range(dates), "1 months"); 
atVector  = seq(r[1], as.POSIXct("2018-01-01 01:00:00 CET"), by = "6 months")
atVector2 = seq(r[1], as.POSIXct("2018-01-01 01:00:00 CET"), by = "12 months")
axis.Date(1, at=atVector, labels=NA)
axis.Date(1, at=atVector2, format="%Y")

grid = as.Date(dates);
polygon(c(grid,rev(grid)),c(unemployment,rep(4,length(grid))),border=NA,col=rgb(t(col2rgb(plotColors[2]))/256,alpha=0.25))
lines(as.Date(dates),unemployment,lwd=1,col=plotColors[2])
abline(v=as.Date("2016-01-01 01:00:00 CET"),lty="dotted")

plot(as.Date(dates),inflation,type="l",ylab="inflation rate", xlab="date",col=plotColors[1],bty="n",xaxt="n", ylim=c(-2,4), xlim=c(as.Date("2015-01-01 01:00:00 CET"),as.Date("2018-01-01 01:00:00 CET")) )

r         = as.POSIXct(range(dates), "1 months"); 
atVector  = seq(r[1], as.POSIXct("2018-01-01 01:00:00 CET"), by = "6 months")
atVector2 = seq(r[1], as.POSIXct("2018-01-01 01:00:00 CET"), by = "12 months")
axis.Date(1, at=atVector, labels=NA)
axis.Date(1, at=atVector2, format="%Y")

grid = as.Date(dates);

upper_inflation = apply(y,2,mean) + 1.96 * apply(y,2,sd)
lower_inflation = apply(y,2,mean) - 1.96 * apply(y,2,sd)

polygon(c(grid[1:13],rev(grid[1:13])),c(inflation[1:13],rep(-2,13)),border=NA,col=rgb(t(col2rgb(plotColors[1]))/256,alpha=0.25))
polygon(c(grid[13:36],rev(grid[13:36])),c(upper_inflation,rev(lower_inflation)),border=NA,col=rgb(t(col2rgb(plotColors[3]))/256,alpha=0.25))
lines(as.Date(dates),inflation,lwd=1,col=plotColors[1])
lines(grid[13:36],apply(y,2,mean),lwd=2,col=plotColors[3])
abline(v=as.Date("2016-01-01 01:00:00 CET"),lty="dotted")

dev.off()
