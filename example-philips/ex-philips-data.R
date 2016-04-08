##############################################################################
##############################################################################
# Example 2.1
# Philips curve modelfor Swedish inflation/unemployment
# Reproduces Figure 2.2
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

d            <- read.table("philipscurve_sweden_1987_2015.csv",sep=",",header=T)
dates        <- as.Date( seq(as.POSIXct("1987-01-01 01:00:00 CET"), as.POSIXct("2015-12-01 01:00:00 CET"), by = "1 months") )
inflation    <- as.numeric( d$inflation    )
unemployment <- as.numeric( d$unemployment )
nMCMC        <- 15000;
burnin       <- 5000;


cairo_pdf("ex-philips-data.pdf",  height = 8, width = 8)

layout(matrix(c(1,2), 2, 1, byrow = TRUE))  
par(mar=c(4,5,1,1))

plot(as.Date(dates),inflation,type="l",ylab="inflation rate", xlab="date",col=plotColors[1],bty="n",xaxt="n", ylim=c(-2,15), xlim=c(as.Date("1987-01-01 01:00:00 CET"),as.Date("2015-12-01 01:00:00 CET")) )

tt <- c(as.Date("1991-01-01 01:00:00 CET"),as.Date("1993-01-01 01:00:00 CET"))
polygon(c(tt,rev(tt)),c(-2,-2,15,15),border=NA,col=rgb(t(col2rgb(plotColors[3]))/256,alpha=0.25))

tt <- c(as.Date("2008-01-01 01:00:00 CET"),as.Date("2010-01-01 01:00:00 CET"))
polygon(c(tt,rev(tt)),c(-2,-2,15,15),border=NA,col=rgb(t(col2rgb(plotColors[3]))/256,alpha=0.25))

r         = as.POSIXct(range(dates), "1 months"); 
atVector  = seq(r[1], as.POSIXct("2016-01-01 01:00:00 CET"), by = "12 months")
atVector2 = seq(r[1], as.POSIXct("2016-01-01 01:00:00 CET"), by = "24 months")
axis.Date(1, at=atVector, labels=NA)
axis.Date(1, at=atVector2, format="%Y")

grid = as.Date(dates);
polygon(c(grid,rev(grid)),c(inflation,rep(-2,length(grid))),border=NA,col=rgb(t(col2rgb(plotColors[1]))/256,alpha=0.25))
lines(as.Date(dates),inflation,lwd=1,col=plotColors[1])


plot(as.Date(dates),unemployment,type="l",ylab="unemployment rate", xlab="date",col=plotColors[2],bty="n",xaxt="n", ylim=c(0,15) , xlim=c(as.Date("1987-01-01 01:00:00 CET"),as.Date("2015-12-01 01:00:00 CET")) )

tt <- c(as.Date("1991-01-01 01:00:00 CET"),as.Date("1993-01-01 01:00:00 CET"))
polygon(c(tt,rev(tt)),c(0,0,15,15),border=NA,col=rgb(t(col2rgb(plotColors[3]))/256,alpha=0.25))

tt <- c(as.Date("2008-01-01 01:00:00 CET"),as.Date("2010-01-01 01:00:00 CET"))
polygon(c(tt,rev(tt)),c(0,0,15,15),border=NA,col=rgb(t(col2rgb(plotColors[3]))/256,alpha=0.25))

r         = as.POSIXct(range(dates), "1 months"); 
atVector  = seq(r[1], as.POSIXct("2016-01-01 01:00:00 CET"), by = "12 months")
atVector2 = seq(r[1], as.POSIXct("2016-01-01 01:00:00 CET"), by = "24 months")
axis.Date(1, at=atVector, labels=NA)
axis.Date(1, at=atVector2, format="%Y")

grid = as.Date(dates);
polygon(c(grid,rev(grid)),c(unemployment,rep(0,length(grid))),border=NA,col=rgb(t(col2rgb(plotColors[2]))/256,alpha=0.25))
lines(as.Date(dates),unemployment,lwd=1,col=plotColors[2])

dev.off()