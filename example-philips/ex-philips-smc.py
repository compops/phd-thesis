##############################################################################
##############################################################################
# Example 3.3
# Particle filtering for Swedish inflation/unemployment
#
# Copyright (c) 2016 Johan Dahlin [ johan.dahlin (at) liu.se ]
# Distributed under the MIT license.
#
##############################################################################
##############################################################################


import numpy   as np
import pandas
from   state   import smc
from   models  import philipscurve_4parameters

##############################################################################
# Arrange the data structures
##############################################################################
sm               = smc.smcSampler();

##############################################################################
# Setup the system
##############################################################################
sys               = philipscurve_4parameters.ssm()
sys.par           = np.zeros((sys.nPar,1))
sys.par[0]        = 0.75520623;
sys.par[1]        = 0.45212950;
sys.par[2]        = 0.01244553;
sys.par[3]        = 0.27753243;
sys.T             = 348;
sys.xo            = 2.0;


##############################################################################
# Load data
##############################################################################

inp = np.loadtxt("data/philipscurve_sweden_1987_2015.csv",delimiter=",",dtype="float",skiprows=1)

kpi = inp[:,0];
aku = inp[:,1];

# Unemployment rate
sys.u = aku;

# Inflation
sys.y = kpi;


##############################################################################
# Setup the parameters
##############################################################################
th               = philipscurve_4parameters.ssm()
th.nParInference = 4;
th.copyData(sys);


##############################################################################
# Setup the SMC algorithm
##############################################################################

sm.genInitialState = True;
sm.nPart           = 100;
sm.resampFactor    = 2.5;
sm.genInitialState = True;
sm.xo              = sys.xo;
th.xo              = sys.xo;

##############################################################################
# Run the particle filter
##############################################################################

sm.bPF(th);
sm.writeToFile();

##############################################################################
# Run the particle filter
##############################################################################

nPartVector = (50, 100, 250, 500, 1000)
llMatrix    = np.zeros(((len(nPartVector),100)))

for ii in range(len(nPartVector)):
    sm.nPart = nPartVector[ii];

    for jj in range(100):
        sm.bPF(th);
        llMatrix[ii,jj] = sm.ll;
        print((ii,jj))

np.savetxt("philips_smc_ll.csv",llMatrix,delimiter=',')

np.mean(llMatrix,axis=1)
np.sqrt(np.var(llMatrix,axis=1))

##############################################################################
# Run the particle filter
##############################################################################

sm.nPart  = 50;
sm.calcGradientFlag = True
sm.filter = sm.bPF;
phiGrid   = np.arange(0.59,1.0,0.01)
llMatrix  = np.zeros(((len(phiGrid),2)))

for ii in range(len(phiGrid)):
    th.par[0] = phiGrid[ii]
    sm.flPS(th);
    llMatrix[ii,:] = ( sm.ll, sm.gradient[0] );
    print((ii))

subplot(2,1,1); plot(phiGrid,llMatrix[:,0])
subplot(2,1,2); plot(phiGrid,llMatrix[:,1])

np.savetxt("philips_smc_llgrid.csv",llMatrix,delimiter=',')

##############################################################################
##############################################################################
# End of file
##############################################################################
##############################################################################
