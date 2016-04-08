##############################################################################
##############################################################################
# Example 1.1
# Particle Metropolis-Hastings to estimate the parameters and state
#
# Copyright (c) 2016 Johan Dahlin [ johan.dahlin (at) liu.se ]
# Distributed under the MIT license.
#
##############################################################################
##############################################################################

import pandas

import numpy            as np
import matplotlib.pylab as plt

from   state   import smc
from   para    import pmh
from   models  import icevarve_4parameters

##############################################################################
# Arrange the data structures
##############################################################################
sm               = smc.smcSampler();
pmh              = pmh.stPMH();

##############################################################################
# Setup the system
##############################################################################
sys               = icevarve_4parameters.ssm()
sys.par           = np.zeros((sys.nPar,1))
sys.par[0]        = 0.953;
sys.par[1]        = 0.034;
sys.par[2]        = 6.420;
sys.par[3]        = 0.285;
sys.T             = 634;
sys.xo            = 0.0;
sys.version       = "standard"
sys.filePrefix    = "icevarve_4parameters";


##############################################################################
# Generate data
##############################################################################
sys.generateData(np.zeros(sys.T),"data/icevarve.txt","y");


##############################################################################
# Setup the parameters
##############################################################################
th               = icevarve_4parameters.ssm()
th.version       = "standard"
th.nParInference = 4;
th.copyData(sys);


##############################################################################
# Run particle filter
##############################################################################
sm.filter          = sm.bPF;
sm.nPart           = 1000;
sm.genInitialState = True;

##############################################################################
# Run PMH
##############################################################################

# Settings for the PMH routine
pmh.initPar       = ( 0.50,  1.00,  5.85462963,  0.20759259 )
pmh.invHessian    = np.matrix([[ 2.635012e-04, -6.279865e-05, -0.002609949, -2.872486e-04 ],
                               [ -6.279865e-05,  3.453545e-05,  0.001444188,  9.898419e-05 ],
                               [ -2.609949e-03,  1.444188e-03,  0.212732355,  9.775651e-03 ],
                               [ -2.872486e-04,  9.898419e-05,  0.009775651,  1.480525e-03 ]]);
pmh.stepSize      = 2.562 / np.sqrt( th.nParInference );
pmh.nIter         = 10000;
pmh.nBurnIn       = 2500;

pmh.dataset        = 0;
pmh.runSampler(sm,sys,th,"pPMH0")
pmh.writeToFile(sm);


foo = sys.par[2] * np.exp( pmh.x[pmh.nBurnIn:pmh.nIter,:] ) / sys.par[3]
np.savetxt("icevarve_pmh_xhats.csv",np.mean(foo,axis=0),delimiter=',')
np.savetxt("icevarve_pmh_xhats_var.csv",np.var(foo,axis=0),delimiter=',')

##############################################################################
##############################################################################
# End of file
##############################################################################
##############################################################################
