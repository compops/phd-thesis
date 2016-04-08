##############################################################################
##############################################################################
# Example 1.1
# Gaussian process regression for ice varve data
#
# Copyright (c) 2016 Johan Dahlin [ johan.dahlin (at) liu.se ]
# Distributed under the MIT license.
#
##############################################################################
##############################################################################

import numpy as np
import pylab as pb
import GPy
import pandas

d = np.loadtxt("data/icevarve.txt")

# Generate some realisations from the GP prior
x       = np.arange(0,634,1)
X       = x[:,None];
Y       = d[:,None];
mu      = np.zeros(634); 

# Fit the GP
k1 = GPy.kern.Bias(1) + GPy.kern.Matern32(1, lengthscale=1);      
m1 = GPy.models.GPRegression(X,Y,k1); 

m1.sum.bias.variance           = 767.366976247
m1.sum.Mat32.lengthscale       = 25
m1.sum.Mat32.variance          = 226.641568122
m1.Gaussian_noise.variance     = 173.380045038

#m1.optimize('bfgs',max_iters=200)
#m1.optimize_restarts(num_restarts = 10, robust=True)

# Evaluate the predictive distribution on a grid
Mup1, var1 = m1.predict( X ); 

# Export data
out = np.hstack((Mup1,Mup1-1.96*np.sqrt(var1),Mup1+1.96*np.sqrt(var1)));
pandas.DataFrame(out).to_csv("ch1-icevarve-posterior.csv");
pandas.DataFrame(X).to_csv("ch1-icevarve-grid.csv");

########################################################################
# End of file
########################################################################
