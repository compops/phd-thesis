##############################################################################
##############################################################################
# Example 3.13
# Particle Metropolis-Hastings for Swedish inflation/unemployment
#
# Copyright (c) 2016 Johan Dahlin [ johan.dahlin (at) liu.se ]
# Distributed under the MIT license.
#
##############################################################################
##############################################################################

import numpy   as np
from   state   import smc
from   para    import pmh
from   models  import philipscurve_4parameters

##############################################################################
# Arrange the data structures
##############################################################################
sm               = smc.smcSampler();
pmh              = pmh.stPMH();

##############################################################################
# Setup the system
##############################################################################
sys               = philipscurve_4parameters.ssm()
sys.par           = np.zeros((sys.nPar,1))
sys.par[0]        = 0.93;
sys.par[1]        = 0.42;
sys.par[2]        = -0.331;
sys.par[3]        = 1.0;
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

sm.filter          = sm.bPF;
sm.genInitialState = True;

sm.nPart           = 100;
sm.resampFactor    = 2.5;
sm.genInitialState = True;
sm.xo              = sys.xo;
th.xo              = sys.xo;

##############################################################################
# Setup the PMH algorithm
##############################################################################

pmh.dataset              = 0;
pmh.nIter                = 15000;
pmh.nBurnIn              = 5000;
pmh.initPar              = ( 0.75732627,  0.44732395,  0.01203858,  0.27690514 );

# Settings for the PMH routine
pmh.stepSize             = 2.562 / np.sqrt(th.nParInference);
pmh.invHessian           = np.matrix([[  7.25303580e-03,  -2.96004273e-03,  -1.85646162e-06,  -5.41850306e-05],
                                      [ -2.96004273e-03,   2.86801751e-02,  -7.38086096e-05,  -5.84915851e-05],
                                      [ -1.85646162e-06,  -7.38086096e-05,   2.80899543e-05,  -6.11368058e-06],
                                      [ -5.41850306e-05,  -5.84915851e-05,  -6.11368058e-06,   1.02259119e-04]])
pmh.nProgressReport      = 100;

##############################################################################
# Run the qPMH2 sampler
##############################################################################

pmh.runSampler(sm, sys, th, "pPMH0");
pmh.writeToFile();

xhatMean = np.mean( pmh.x[ pmh.nBurnIn:pmh.nIter, : ] , axis=0 )
xhatSTD  = np.sqrt(np.var( pmh.x[ pmh.nBurnIn:pmh.nIter, : ] , axis=0 ))

np.savetxt("philips_nairu_mean.csv",xhatMean,delimiter=',')
np.savetxt("philips_nairu_sd.csv",xhatSTD,delimiter=',')

##############################################################################
##############################################################################
# End of file
##############################################################################
##############################################################################
