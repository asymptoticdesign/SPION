#!/usr/bin/python

import scipy
import numpy
import pylab
import spionlib


def example1d():
    None

def example2d():
    None

def hysteresis_example():
    None

spions = []

spions.append(('Fe',spionlib.Fe_SPION(diameter = 10.0e-9)))
spions.append(('Mg',spionlib.Mg_SPION(diameter = 10.0e-9)))
spions.append(('Co',spionlib.Co_SPION(diameter = 10.0e-9)))
spions.append(('Ni',spionlib.Ni_SPION(diameter = 10.0e-9)))
spions.append(('Mn',spionlib.Mn_SPION(diameter = 10.0e-9)))
spions.append(('Zn',spionlib.Zn_SPION(diameter = 10.0e-9)))

for spion in spions:
    tau = spion[1].computeNeel()
    print spion[1].computeResonance(tau)/1.0e6,"MHz"

H,M = spions[0][1].computeHysteresis(1.0e-2,1.0e3)

pylab.plot(H,M)
pylab.show()
