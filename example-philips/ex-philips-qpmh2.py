##############################################################################
##############################################################################
# Example 4.1
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
th.nParInference = 3;
th.copyData(sys);


##############################################################################
# Setup the SMC algorithm
##############################################################################

sm.filter          = sm.bPF;
sm.smoother        = sm.flPS;
sm.genInitialState = True;

sm.nPart           = 100;
sm.resampFactor    = 2.5;
sm.genInitialState = True;
sm.fixedLag        = 12;
sm.xo              = sys.xo;
th.xo              = sys.xo;

##############################################################################
# Setup the PMH algorithm
##############################################################################

pmh.dataset              = 0;
pmh.nIter                = 10000;
pmh.nBurnIn              = 2500;
pmh.initPar              = ( 0.75732627,  0.44732395,  0.01203858,  0.27690514 );

# Settings for the qPMH2 routine

pmh.stepSize             = 1.0;
pmh.epsilon              = 10000;
pmh.memoryLength         = 100;
pmh.nProgressReport      = 100;

# Use the hybrid form of the sampler (replacing a non-PSD Hessian with the
# estimated posterior covariance using the last iterations of the burn-in)
pmh.makeHessianPSDmethod          = "hybrid";
pmh.PSDmethodhybridSamps          = 1500;
pmh.qPMHadaptInitialHessian       = True;


##############################################################################
# Run the qPMH2 sampler
##############################################################################

pmh.runSampler(sm, sys, th, "qPMH2");

np.savetxt("philips_post_qpmh2.csv",pmh.th,delimiter=',')
np.savetxt("philips_post_qpmh2zv.csv",pmh.thzv,delimiter=',')

##############################################################################
##############################################################################
# End of file
##############################################################################
##############################################################################
