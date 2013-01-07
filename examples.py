#!/usr/bin/python

import scipy
import numpy
import pylab
import spionlib


def example1d():
    """
    Computes the specific loss power considering only neel relaxation (slpN), only brown relaxation (slpB), and all relaxation mechanisms (slpT).
    """
    #Add field parameters -- 1 MHz and 1000 A/m
    spion[1].coupleField(1.0e6,1.0e4)
    d = scipy.linspace(1.0e-9,50.0e-9,100)
    slpT = spion[1].computeSLP(d=d)
    slpB = spion[1].computeSLP(tau=spion[1].computeBrown(d=d),d = d)
    slpN = spion[1].computeSLP(tau=spion[1].computeNeel(d=d),d = d)

    pylab.plot(d/1.0,slpT,color='black')
    pylab.plot(d/1.0,slpB,color='red')
    pylab.plot(d/1.0,slpN,color='blue')
    pylab.show()

def example2d():
    None

def hysteresis_example():
    None

spions = []

spions.append(('Fe',spionlib.Fe3O4_SPION(diameter = 10.0e-9)))
spions.append(('Mg',spionlib.Mg_SPION(diameter = 10.0e-9)))
spions.append(('Co',spionlib.Co_SPION(diameter = 10.0e-9)))
spions.append(('Ni',spionlib.Ni_SPION(diameter = 10.0e-9)))
spions.append(('Mn',spionlib.Mn_SPION(diameter = 10.0e-9)))
spions.append(('Zn',spionlib.Zn_SPION(diameter = 10.0e-9)))

for spion in spions:
    spion[1].coupleCarrier(0.001)
    spion[1].computeResonance()
    print spion[1].resonance/1.0e6

#H,M = spions[0][1].computeHysteresis(1.0e-2,1.0e-10)
#pylab.plot(H,M)
#pylab.show()

example1d()


