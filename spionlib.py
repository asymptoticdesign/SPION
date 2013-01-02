"""
Title: SPIONlib
Description: An object-oriented library for computing parameters of superparamagnetic iron oxide nanoparticles.
"""

import scipy
import numpy

#physical constants
u0 = scipy.pi*4.0e-7 #H / m
kb = 1.38e-23 #joules / kelvin
T0 = 298.0 #kelvin
kbT = 4.114e-21 #in joules

eta_h2o = 0.001
eta_toluene = 0.00059
rho_fe = 5.175e6 #g / m^3

class SPION:
    """
    A class that represents a SPION.  The idea of the class is to give the nanoparticle a set of 'default' values, and then every computation can optionally accept a vector argument via scipy.  Any non-vector arguments are taken from the default values.
    """

    def __init__(self, diameter, anisotropy, magnetization, density=5.175e6):
        """
        Input parameters to initialize a SPION of a particular material.
        diameter: in nanometers
        anisotropy constant: in J / m^3
        specific magnetization (sigma): A m^2 / g
        density = 5.175e6 g / m^3
        """
        #check if all inputs are floats (and not lists)
        for i in [diameter, anisotropy, magnetization]:
            if !isinstance.(i, float):
                raise TypeError('Invalid input.  All initial inputs must be floats.')
            
        self.diameter = diameter
        self.volume = 4/3 * scipy.pi * (self.diameter/2)**3
        self.K = anisotropy
        self.sigma = magnetization
        self.density = density

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
            raise Exception('Viscosity (eta) not set.  Either pass viscosity value or set a default using coupleCarrier()')
        
        return (eta * V) / kbT

    def computeRelaxation(self, eta = self.eta, V = self.volume, K = self.K):
        """
        Compute the relaxation time (in seconds) for the particle in a given carrier fluid.
        """
        if(!eta):
            raise Exception('Viscosity (eta) not set.  Either pass viscosity value or set a default using coupleCarrier()')
        
        neel = self.computeNeel(K,V)
        brown = self.computeBrown(eta,V)
        return neel*brown / (neel + brown)

    def computeResonance(self, tau):
        """
        Computes the 'relaxation resonance' frequency (in Hz) at which the most magnetic losses occur.
        """
        return 1.0 / (2*pi*tau)

    def computeResonance(self, K = self.K, V = self.volume, eta = self.eta):
        """
        Computes the 'relaxation resonance' frequency (in Hz) at which the most magnetic losses occur.
        """
        if(!eta):
            raise Exception('Viscosity (eta) not set.  Either pass viscosity value or set a default using coupleCarrier()')
        
        return 1.0 / (2*pi*self.computeRelaxation(eta,V,K))

    def __prefactor__(self,fieldAmplitude, fieldFrequency, sigma):
        prefactor = u0 * scipy.pi * sigma * fieldAmplitude * fieldFrequency
        return prefactor

    def __xi__(self,fieldAmplitude, sigma, volume):
        xi = (u0 * sigma * self.density * fieldAmplitude * volume)/ kbT
        return xi

    def __quadSuscept__(self, tau, fieldFreq):
        quadSuscept = (2*scipy.pi*fieldFreq*tau) / (1 + (2*scipy.pi*fieldFreq*tau)**2)
        return quadSuscept

    def computeSLP(self, tau, K = self.K, V = self.volume, sigma = self.sigma, fieldAmp = self.fieldAmp, fieldFreq = self.fieldFreq):
        """
        Compute the specific loss power in W / g*s.  'tau' is a relaxation time constant as output by self.computeNeel, self.computeBrown, or self.computeRelaxation.  
        """
        prefactor = self.__prefactor__(fieldAmp, fieldFreq, sigma)
        xi = self.__xi__(fieldAmp, sigma, volume)
        quadSuscept = self.__quadSuscept__(tau,fieldFreq)
        slp = prefactor * quadSuscept * (scipy.tanh(xi)**-1 - (1/xi))
        return slp
        
#Now various other functions
def make2dVectors(vector1, vector2):
    """
    create a meshgrid using the two input vectors for two-parameter simulation.
    """
    #check that the vector are correct length, etc.
    if (len(vector1) != len(vector2)):
        raise Exception("Both input vectors must be the same length.")
    
    output1, output2 = numpy.meshgrid(vector1, vector2)
    return output1, output2

#Predefined SPIONs
def Fe_SPION(diameter=10.0e-9):
    return SPION(diameter, 1.4e4, 0.092)

def Co_SPION(diameter=10.0e-9):
    return SPION(diameter, 1.8e4, 0.08)

def Mn_SPION(diameter=10.0e-9):
    return SPION(diameter, 3.3e3, 0.08)

def Ni_SPION(diameter=10.0e-9):
    return SPION(diameter, 3.3e3, 0.05)

def Mg_SPION(diameter=10.0e-9):
    return SPION(diameter, 2.0e4, 0.03)

def Zn_SPION(diameter=10.0e-9):
    return SPION(diameter, 4.6e3, 0.065)
