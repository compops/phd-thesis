##############################################################################
##############################################################################
# Helpers for distributions and the gradients and Hessian of the log-pdfs
#
# Copyright (c) 2016 Johan Dahlin [ johan.dahlin (at) liu.se ]
# Distributed under the MIT license.
#
##############################################################################
##############################################################################

import numpy                 as     np
import scipy                 as     sp
import math

##############################################################################
# Half Cauchy pdf, gradient and Hessian of log-pdf
##############################################################################

def halfCauchyPDF(x, x0, gamma):
    return 2.0 / ( np.pi * np.sqrt( gamma ) ) * ( 1.0 + ( x - x0 )**2 / gamma );

def halfCauchyLogPDF(x, x0, gamma):
    return np.log( 2.0 ) - np.log( 1.0 + ( x - x0 )**2 / gamma ) - np.log( np.pi * np.sqrt( gamma ) );

##############################################################################
# Gamma pdf, gradient and Hessian of log-pdf
##############################################################################

def gammaPDF(x, a, b):
    # Shape and rate parametrisation
    return b**a / sp.special.gamma(a) * x**(a-1.0) * np.exp(-b*x)

def gammaLogPDF(x, a, b):
    # Shape and rate parametrisation
    return a * np.log(b) + (a-1.0) * np.log(x) - b * x - sp.special.gammaln(a)

def gammaLogPDFgradient (x, a, b):
    # Shape and rate parametrisation
    return ( a - 1.0 ) / x - b;

def gammaLogPDFhessian (x, a, b):
    # Shape and rate parametrisation
    return - ( a - 1.0 ) / ( x**2 );

##############################################################################
# Inverse Gamma pdf, gradient and Hessian of log-pdf
##############################################################################

def invGammaPDF(x, a, b):
    # Shape and rate parametrisation
    return b**a / sp.special.gamma(a) * x**(-a-1.0) * np.exp(-b/x)

def invGammaLogPDF(x, a, b):
    # Shape and rate parametrisation
    return a * np.log(b) + (-a-1.0) * np.log(x) - (b / x) - sp.special.gammaln(a)

def invGammaLogPDFgradient (x, a, b):
    # Shape and rate parametrisation
    return ( - a - 1.0 ) / x + b / ( x**2 );

def invGammaLogPDFhessian (x, a, b):
    # Shape and rate parametrisation
    return - ( - a - 1.0 ) / ( x**2 ) - 2.0 * b / ( x**3 );

##############################################################################
# Beta pdf, gradient and Hessian of log-pdf
##############################################################################

def betaPDF(x, a, b):
    return sp.special.gamma(a+b) / ( sp.special.gamma(a) + sp.special.gamma(b) ) * x**(a-1.0) * (1.0-x)**(b-1.0);

def betaLogPDF(x, a, b):
    return sp.special.gammaln(a+b) - np.log( sp.special.gamma(a) + sp.special.gamma(b) ) + (a-1.0)*np.log(x) * (b-1.0) * np.log(1.0-x);

def betaLogPDFgradient (x, a, b):
    return ( a - 1.0 ) / x + ( 1.0 - b ) / ( 1.0 - x );

def betaLogPDFhessian (x, a, b):
    return - ( a - 1.0 ) / ( x**2 ) + ( 1.0 - b ) / ( ( 1.0 - x )**2 );

##############################################################################
# Normal pdf, gradient and Hessian of log-pdf
##############################################################################

def normalPDF(x, a, b):
    return 1.0 / np.sqrt( 2 * np.pi * b**2 ) * np.exp( -0.5/(b**2) * (x-a)**2 );

def normalLogPDF(x, a, b):
    return -0.5 * np.log( 2 * np.pi * b**2 ) -0.5/(b**2) * (x-a)**2;

def normalLogPDFgradient (x, a, b):
    return (a - x)/b**2;

def normalLogPDFhessian (x, a, b):
    return -1.0/b**2

##############################################################################
# Multivariate normal pdf
##############################################################################

def MultivariateNormalPDF(x, mu, S):
    nx = len(S)
    norm_coeff = nx * np.log( 2.0 * np.pi ) + np.linalg.slogdet(S)[1]
    err = x-mu

    numerator = np.dot( np.dot(err,np.linalg.pinv(S)),err.transpose())
    return np.exp( -0.5*(norm_coeff+numerator) )

def MultivariateNormalLogPDF(x, mu, S):
    nx = len(S)
    norm_coeff = nx * np.log( 2.0 * np.pi ) + np.linalg.slogdet(S)[1]
    err = x-mu

    numerator = np.dot( np.dot(err,np.linalg.pinv(S)),err.transpose())
    return -0.5*(norm_coeff+numerator)


##############################################################################
##############################################################################
# End of file
##############################################################################
##############################################################################
