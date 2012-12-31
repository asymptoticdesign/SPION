"""
Title: SPIONlib
Description: An object-oriented library for computing parameters of superparamagnetic iron oxide nanoparticles.
"""

import scipy

#physical constants
pi = scipy.pi
u0 = pi*4.0e-7 #H / m
kb = 1.38e-23 #joules / kelvin
T0 = 298.0 #kelvin
kbT = 4.114e-21 #in joules

eta_h2o = 0.001
eta_t = 0.00059
rho_fe = 5.175e6 #g / m^3

class SPION:
    """
    A class that represents a SPION.  The idea of the class is to give the nanoparticle a set of 'default' values, and then every computation can optionally accept a vector argument via scipy.  Any non-vector arguments are taken from the default values.
    """

    def __init__(self, diameter, anisotropy, magnetization):
        """
        Input parameters to initialize a SPION of a particular material.
        diameter: in nanometers
        anisotropy constant: in J / m^3
        specific magnetization (sigma): A m^2 / g
        """
        #check if all inputs are floats (and not lists)
        for i in [diameter, anisotropy, magnetization]:
            if !isinstance.(i, float):
                raise TypeError('Invalid input.  All initial inputs must be floats.')
            
        self.diameter = diameter
        self.volume = 4/3 * scipy.pi * (self.diameter/2)**3
        self.K = anisotropy
        self.sigma = magnetization

    def coupleField(self, fieldFreq, fieldAmp):
        """
        Add default field parameters for future computations.  Field frequency should be in Hertz and field amplitude should be in Amps / meter
        """
        for i in [fieldFreq, fieldAmp]:
            if !isinstance.(i, float):
                raise TypeError('Invalid input.  All field parameters must be floats.')
        self.fieldFreq = fieldFreq
        self.fieldAmp = fieldAmp

    def coupleCarrier(self, eta):
        """
        Add a default carrier fluid.
        """
        if !isinstance.(eta, float):
            raise TypeError('Invalid input.  Viscosity (eta) must be a float.')
        self.eta = eta

    def computeNeel(self, K = self.K, V = self.volume):
        """
        Compute the Neel Relaxation time constant in seconds.
        """
        return 1.0e-9 * scipy.exp(K*V/(kbT))

    def computeBrown(self, eta = self.eta, V = self.volume):
        """
        Compute the Brown Relaxation time constant in seconds.
        """
        if(!eta):
            raise TypeError('Viscosity (eta) not set.  Either pass viscosity value or set a default using coupleCarrier()')
        
        return (eta * V) / kbT

    def computeRelaxation(self, eta = self.eta, V = self.volume, K = self.K):
        """
        Compute the relaxation time (in seconds) for the particle in a given carrier fluid.
        """
        if(!eta):
            raise TypeError('Viscosity (eta) not set.  Either pass viscosity value or set a default using coupleCarrier()')
        
        neel = self.computeNeel(K,V)
        brown = self.computeBrown(eta,V)
        return neel*brown / (neel + brown)    

    def computeResonance(self, K = self.K, V = self.volume, eta = self.eta):
        """
        Computes the 'relaxation resonance' frequency (in Hz) at which the most magnetic losses occur.
        """
        if(!eta):
            raise TypeError('Viscosity (eta) not set.  Either pass viscosity value or set a default using coupleCarrier()')
        
        return 1.0 / (2*pi*self.computeRelaxation(eta,V,K))

    def computeSLP(self, tau, K = self.K, V = self.volume, sigma = self.sigma):
        """
        Compute the specific loss power in W / g*s.  'tau' is a relaxation time constant as output by self.computeNeel, self.computeBrown, or self.computeRelaxation.
        """
        
#Now various other functions
