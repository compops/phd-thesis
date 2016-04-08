##############################################################################
##############################################################################
# Example 3.13
# Particle filtering for Swedish inflation/unemployment
# Reproduces Figure 3.6
#
# Copyright (c) 2016 Johan Dahlin [ johan.dahlin (at) liu.se ]
# Distributed under the MIT license.
#
##############################################################################
##############################################################################


# Setup plot colors
library("RColorBrewer")
plotColors = brewer.pal(7, "Dark2");

d            <- read.table("philipscurve_sweden_1987_2015.csv",sep=",",header=T)
dates        <- as.Date( seq(as.POSIXct("1987-01-01 01:00:00 CET"), as.POSIXct("2015-12-01 01:00:00 CET"), by = "1 months") )
inflation    <- as.numeric( d$inflation    )
unemployment <- as.numeric( d$unemployment )

dpmh         <- read.table("pmh0/philips_th_posterior.csv",header=T,sep=",");
xhats        <- as.numeric( read.table("pmh0/philips_nairu_mean.csv",header=F,sep=",")$V1 );
xhatm        <- as.numeric( read.table("smc/state_bPF_N100.csv",header=T,sep=",")$xhatf );

cairo_pdf("ex-philips-state.pdf",  height = 6, width = 8)

layout(matrix(1:4, 2, 2, byrow = TRUE))  
par(mar=c(4,5,1,1))

plot(as.Date(dates),xhatm,type="n",ylab="NAIRU (bPF)", xlab="date",col=plotColors[2],bty="n",xaxt="n", ylim=c(1,3) , xlim=c(as.Date("1987-01-01 01:00:00 CET"),as.Date("2015-12-01 01:00:00 CET")) )

r         = as.POSIXct(range(dates), "1 months"); 
atVector  = seq(r[1], as.POSIXct("2016-01-01 01:00:00 CET"), by = "24 months")
atVector2 = seq(r[1], as.POSIXct("2016-01-01 01:00:00 CET"), by = "48 months")
axis.Date(1, at=atVector, labels=NA)
axis.Date(1, at=atVector2, format="%Y")

grid = as.Date(dates);
polygon(c(grid,rev(grid)),c(xhatm,rep(1,length(grid))),border=NA,col=rgb(t(col2rgb(plotColors[3]))/256,alpha=0.25))
lines(as.Date(dates)[-1],xhatm[-1],lwd=1,col=plotColors[3])

plot(as.Date(dates),unemployment-xhatm,type="n",ylab="unemployment gap (bPF)", xlab="date",col=plotColors[4],bty="n",xaxt="n", ylim=c(-2,12) , xlim=c(as.Date("1987-01-01 01:00:00 CET"),as.Date("2015-12-01 01:00:00 CET")) )

r         = as.POSIXct(range(dates), "1 months"); 
atVector  = seq(r[1], as.POSIXct("2016-01-01 01:00:00 CET"), by = "24 months")
atVector2 = seq(r[1], as.POSIXct("2016-01-01 01:00:00 CET"), by = "48 months")
axis.Date(1, at=atVector, labels=NA)
axis.Date(1, at=atVector2, format="%Y")

polygon(c(grid,rev(grid)),c(unemployment-xhatm,rep(-2,length(grid))),border=NA,col=rgb(t(col2rgb(plotColors[4]))/256,alpha=0.25))
lines(as.Date(dates),unemployment-xhatm,lwd=1,col=plotColors[4])



plot(as.Date(dates),xhats,type="l",ylab="unemployment rate (PMH)", xlab="date",col=plotColors[2],bty="n",xaxt="n", ylim=c(1,3) , xlim=c(as.Date("1987-01-01 01:00:00 CET"),as.Date("2015-12-01 01:00:00 CET")) )

r         = as.POSIXct(range(dates), "1 months"); 
atVector  = seq(r[1], as.POSIXct("2016-01-01 01:00:00 CET"), by = "24 months")
atVector2 = seq(r[1], as.POSIXct("2016-01-01 01:00:00 CET"), by = "48 months")
axis.Date(1, at=atVector, labels=NA)
axis.Date(1, at=atVector2, format="%Y")

polygon(c(grid,rev(grid)),c(xhats,rep(1,length(grid))),border=NA,col=rgb(t(col2rgb(plotColors[5]))/256,alpha=0.25))
lines(as.Date(dates),xhats,lwd=1,col=plotColors[5])

plot(as.Date(dates),unemployment-xhats,type="l",ylab="unemployment gap (PMH)", xlab="date",col=plotColors[7],bty="n",xaxt="n", ylim=c(-2,12) , xlim=c(as.Date("1987-01-01 01:00:00 CET"),as.Date("2015-12-01 01:00:00 CET")) )

r         = as.POSIXct(range(dates), "1 months"); 
atVector  = seq(r[1], as.POSIXct("2016-01-01 01:00:00 CET"), by = "24 months")
atVector2 = seq(r[1], as.POSIXct("2016-01-01 01:00:00 CET"), by = "48 months")
axis.Date(1, at=atVector, labels=NA)
axis.Date(1, at=atVector2, format="%Y")

polygon(c(grid,rev(grid)),c(unemployment-xhats,rep(-2,length(grid))),border=NA,col=rgb(t(col2rgb(plotColors[7]))/256,alpha=0.25))
lines(as.Date(dates),unemployment-xhats,lwd=1,col=plotColors[7])

dev.off()