##############################################################################
##############################################################################
# Example 2.2
# US Supreme court ideological leaning
# Reproduces Figure 2.3
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

for (ii in 1:9){
  d[,ii] = d[,ii] * ii
}

cairo_pdf("ex-supreme-data.pdf", height = 8, width = 8)
layout(matrix(1, 1, 1, byrow = TRUE))  
par(mar=c(4,5,1,1))
image(1:nCases,1:9,d,xlab="case",ylab="",
      col=rgb(t(col2rgb(plotColors2))/256,alpha=0.50),
      yaxt="n",xaxt="n",bty='n')
axis(1,at=c(1,nCases),cex.axis=0.75)
axis(2,at=seq(1,9,1),labels=colnames(d),las=1,cex.axis=0.75,lwd=0)
dev.off() 
