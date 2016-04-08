##############################################################################
##############################################################################
# Model specification
# Simple Philips curve model with rational inflation expectations
#
# Copyright (c) 2016 Johan Dahlin [ johan.dahlin (at) liu.se ]
# Distributed under the MIT license.
#
##############################################################################
##############################################################################

#=============================================================================
# Model structure
#=============================================================================
#
# untt = par[1] * unt + par[2] * ( 1 + exp(-unt))**(-1) + 1.0/(1.0+exp(abs(ut - unt))) * vt
# ptt  = pt + par[3] * ( utt - unt ) + par[4] * et
# ( vt, et ) ~iid N(0,1)

import numpy          as     np
from   scipy.stats    import norm
from   models_helpers import *
from   models_dists   import *


class ssm(object):
    #=========================================================================
    # Define model settings
    #=========================================================================
    nPar         = 4;
    par          = np.zeros(nPar);
    modelName    = "Philips curve with rational expectations"
    filePrefix   = "philips-rational";
    supportsFA   = False;
    nParInfernce = None;
    nQInference  = None;
    version       = "standard"

    #=========================================================================
    # Define Jacobian and parameter transforms
    #=========================================================================
    def Jacobian( self ):
        if (self.version == "tanhexp"):
            if ( self.nParInference == 1 ):
                return np.log( 1.0 - self.par[0]**2 )
            if ( self.nParInference == 4 ):
                return np.log( 1.0 - self.par[0]**2 ) + np.log( self.par[3] )
        else:
            return 0.0;

    def transform(self):
        if (self.version == "tanhexp"):
            self.par[0] = np.tanh( self.par[0] );

            if ( self.nParInference > 3 ):
                self.par[3] = np.exp ( self.par[3] );
        return None;

    def invTransform(self):
        if (self.version == "tanhexp"):
            self.par[0] = np.arctanh( self.par[0] );

            if ( self.nParInference > 3 ):
                self.par[3] = np.log    ( self.par[3] );
        return None;

    #=========================================================================
    # Define the model
    #=========================================================================
    def generateInitialState( self, nPart ):
        return np.abs( 2.0 + 2.0 * np.random.normal(size=(1,nPart)) );

    def generateState(self, xt, tt):

        # Generate particles
        part1       = self.par[0] * xt;
        part2       = self.par[1] / ( 1.0 + np.exp(-xt)  );
        part3       = 1.0 / ( 1.0 + np.exp( -np.abs( self.u[tt] - xt ) ) );
        return np.abs( part1 + part2 + part3 * np.random.randn(1,len(xt)) );

    def evaluateState(self, xtt, xt, tt):
        part1       = self.par[0] * xt;
        part2       = self.par[1] / ( 1.0 + np.exp(-xt)  );
        part3       = 1.0 / ( 1.0 + np.exp( -np.abs( self.u[tt] - xt ) ) );
        return norm.pdf( xtt, part1 + part2, part3 );

    def generateObservation(self, xt, tt):
        if ( tt == 0 ):
            return self.par[2] * ( self.u[tt] - xt ) + self.par[3] * np.random.randn(1,len(xt));
        else:
            return self.y[tt-1] + self.par[2] * ( self.u[tt] - xt ) + self.par[3] * np.random.randn(1,len(xt));

    def evaluateObservation(self, xt, tt):
        if ( tt == 0 ):
            return 1.0;
        else:
            return norm.logpdf(self.y[tt], self.y[tt-1] + self.par[2] * ( self.u[tt] - xt ), self.par[3] );

    def generateInitialStateRV( self, nPart, u ):
        return np.abs( 2.0 + 2.0 * u[ range(1,nPart+1) ] )

    def generateStateRV(self, xt, tt, u):
        # Generate particles
        part1       = self.par[0] * xt;
        part2       = self.par[1] / ( 1.0 + np.exp(-xt)  );
        part3       = 1.0 / ( 1.0 + np.exp( -np.abs( self.u[tt] - xt ) ) );
        return np.abs( part1 + part2 + part3 * u[ range(1,len(xt)+1) ] );

    #=========================================================================
    # Define gradients of logarithm of complete data-likelihood
    #=========================================================================
    def Dparm(self, xtt, xt, pu, pw, tt):
        nOut = len(xtt);
        gradient = np.zeros(( nOut, self.nParInference ));
        px    = xtt - self.par[0] * xt - self.par[1] * ( 1.0 + np.exp( -xt ) )**(-1);
        Q2    = ( 1.0 / ( 1.0 + np.exp( -np.abs( self.u[tt] - xt ) ) ) )**(-2);
        py    = self.y[tt] - self.y[tt-1] - self.par[2] * ( self.u[tt] - xt );
        R1    = self.par[3]**(-1);
        R2    = self.par[3]**(-2);
        R3    = self.par[3]**(-3);

        if ( self.version == "standard" ):
            for v1 in range(0,self.nParInference):
                if v1 == 0:
                    gradient[:,v1] = Q2 * px * xt;
                elif v1 == 1:
                    gradient[:,v1] = Q2 * px * ( 1.0 + np.exp( -xt ) )**(-1);
                elif v1 == 2:
                    gradient[:,v1] = R2 * py * ( self.u[tt] - xt );
                elif v1 == 3:
                    gradient[:,v1] = R3 * py**2 - R1;
                else:
                    gradient[:,v1] = 0.0;
        else:
            for v1 in range(0,self.nParInference):
                if v1 == 0:
                    gradient[:,v1] = Q2 * px * xt * ( 1.0 - self.par[0]**2 );
                elif v1 == 1:
                    gradient[:,v1] = Q2 * px * ( 1.0 + np.exp( -xt ) )**(-1);
                elif v1 == 2:
                    gradient[:,v1] = R2 * py * ( self.u[tt] - xt );
                elif v1 == 3:
                    gradient[:,v1] = R2 * py**2 - 1.0;
                else:
                    gradient[:,v1] = 0.0;

        return(gradient);

    #=========================================================================
    # Define Hessians of logarithm of complete data-likelihood
    #=========================================================================
    def DDparm(self, xtt, xt, st, at, tt):

        nOut = len(xtt);
        hessian = np.zeros( (nOut, self.nParInference,self.nParInference) );
        return(hessian);

    #=========================================================================
    # Define hard priors for the PMH sampler
    #=========================================================================
    def priorUniform(self):
        out = 1.0;

        if( np.abs( self.par[0] ) > 1.0 ):
            out = 0.0;

        if( self.par[3] < 0.0 ):
            out = 0.0;

        return( out );

    #=========================================================================
    # Define log-priors for the PMH sampler
    #=========================================================================
    def prior(self):
        out = 0.0;

        # Truncated normal prior for phi (truncation by hard prior)
        if ( self.nParInference >= 1 ):
            out += normalLogPDF( self.par[0], 0.8, 0.1 );

        # Normal prior for alpha
        if ( self.nParInference >= 2 ):
            out += normalLogPDF( self.par[1], 0.5, 0.2 );

        # Normal prior for beta
        if ( self.nParInference >= 3 ):
            out += normalLogPDF( self.par[2], 0.0, 0.1 );

        # Gamma prior for sigma
        if ( self.nParInference >= 4 ):
            out += gammaLogPDF( self.par[3], a=2.0, b=4.0 );

        return out;

    #=========================================================================
    # Define gradients of log-priors for the PMH sampler
    #=========================================================================
    def dprior1(self,v1):

        if ( v1 == 0 ):
            # Truncated normal prior for phi (truncation by hard prior)
            return normalLogPDFgradient( self.par[0], 0.8, 0.1 );
        elif ( v1 == 1):
            # Normal prior for alph
            return normalLogPDFgradient( self.par[1], 0.5, 0.2 );
        elif ( v1 == 2):
            # Normal prior for beta
            return normalLogPDFgradient( self.par[2], 0.0, 0.1 );
        elif ( v1 == 3):
            # Gamma prior for sigma
            return gammaLogPDFgradient( self.par[3], a=2.0, b=4.0 );
        else:
            return 0.0;

    #=========================================================================
    # Define hessians of log-priors for the PMH sampler
    #=========================================================================
    def ddprior1(self,v1,v2):

        if ( ( v1 == 0 ) & ( v1 == 0 ) ):
            # Truncated normal prior for phi (truncation by hard prior)
            return normalLogPDFhessian( self.par[0], 0.8, 0.1 );
        elif ( ( v1 == 1 ) & ( v1 == 1 ) ):
            # Normal prior for alpha
            return normalLogPDFhessian( self.par[1], 0.5, 0.2 );
        elif ( ( v1 == 2 ) & ( v1 == 2 ) ):
            # Normal prior for beta
            return normalLogPDFhessian( self.par[2], 0.0, 0.1 );
        elif ( ( v1 == 3 ) & ( v1 == 3 ) ):
            # Gamma prior for sigma
            return gammaLogPDFhessian( self.par[3], a=2.0, b=4.0 );
        else:
            return 0.0;

    #=========================================================================
    # Sample priors
    #=========================================================================
    def samplePrior(self):

        out = np.zeros( self.nParInference )

        # Truncated normal prior for phi (truncation by hard prior)
        if ( self.nParInference >= 1 ):
            uu = 1.2;
            while ( np.abs(uu) > 1.0):
                uu = np.random.normal( 0.8, 0.1 );
            out[0] = uu;

        # Normal prior for alpha
        if ( self.nParInference >= 2 ):
            out[1] = np.random.normal( 0.5, 0.2 );

        # Normal prior for beta
        if ( self.nParInference >= 3 ):
            out[3] = np.random.normal( 0.0, 0.1 );

        # Gamma prior for sigma
        if ( self.nParInference >= 4 ):
            out[4] = np.random.gamma( shape=2.0, scale=4.0 );

        return out;

    #=========================================================================
    # Define standard methods for the model struct
    #=========================================================================

    # Standard operations on struct
    copyData                = template_copyData;
    storeParameters         = template_storeParameters;
    returnParameters        = template_returnParameters

    # Standard data generation for this model
    generateData            = template_generateData;

    # No faPF available for this model
    generateStateFA         = empty_generateStateFA;
    evaluateObservationFA   = empty_evaluateObservationFA;
    generateObservationFA   = empty_generateObservationFA;

    # No EM algorithm available for this model
    Qfunc                   = empty_Qfunc;
    Mstep                   = empty_Mstep;
