import math
import numpy as np

from scipy.integrate import simps
from scipy.integrate import romberg

from spec_tools.spec_gen import fcoul as fc


#All fundamental constant values taken from NIST website

def spectrum_data_complete(gamma, maxGamma, bFinput, bWMinput=1.41e-3):
    """
    Calculates non-normalized beta decay probability
    for input values of the Lorentz factor gamma, the corresponding
    value of gamma at the high-energy endpoint of the spectrum,
    and the Fierz interference coefficient b.
    The weak magnetism coefficient bWM is an optional input parameter.
    """
    data = (fc.coulomb_function(3., gamma, 6.) * (maxGamma-gamma)**2 * gamma
            * pow(gamma*gamma-1.0,0.5) *
            (1.0 + bFinput / gamma + bWMinput * (2.0 * gamma - maxGamma -1.0 /gamma)))

    return data


def spectrum_probability(energy, b=0):
    """
    Calculates non-normalized beta decay probability for input value
    of electron kinetic energy in eV, and Fierz coefficient b,
    assuming the spectrum endpoint value is that of He6 (K=3.508e6 eV)
    """
    
    me = 5.10998950e5 #Electron rest mass, in eV.
    gamma = energy / me + 1 #convert from kinetic energy in eV to gamma
    maxGamma = 3.508e6 / me + 1 #He6 endpoint gamma
    prob = 0

    if (energy > 3.508e6 or energy < 0.0024):
        print("WARNING: Energy {} out of range".format(energy))
        print("Choose energy between 0.0024 - 3.508e6 eV")
        print()
        return 0

    prob = spectrum_data_complete(gamma, maxGamma, b)
    return prob


def spectrum_prob_norm(b=0):
    """
    Calculates integrated, non-normalized beta decay probability of
    He6 for the full range of He6 decay energies. Used for
    normalization purposes.
    """

    me = 5.10998950e5 #Electron rest mass, in eV.
    maxGamma = 3.508e6/me + 1
    minGamma = 0.0024/me+1

    maxRange =  maxGamma - minGamma

    def function(gamma):
        return spectrum_data_complete(gamma,maxGamma,b)

    normalization = romberg(function,minGamma,maxGamma,divmax=20)

    return normalization


def spectrum_prob_range(min_energy, max_energy, b=0, normalization=0):
    """
    Calculates integrated, non-normalized beta decay probability of He6
    for input range of energies in eV.
    """
    
    me = 5.10998950e5 #Electron rest mass, in eV.
    endGamma = 3.508e6/me + 1 #He6 endpoint gamma
    minGamma = min_energy/me + 1
    maxGamma = max_energy/me + 1
    maxRange=(3.508e6/me + 1)-(0.0024/me + 1)

    if (maxGamma < minGamma):
        print("Maximum Lorentz gamma must be greater than minimum")
        return 0

    if (max_energy > 3.508e6 or min_energy < 0.0024):
        print("Energies out of range")
        print("Min Energy:",min_energy,"eV, Max Energy:",max_energy,"eV")
        print("Choose energies between 0.0024 - 3.508e6 eV")
        return 0

    def function(gamma):
        return spectrum_data_complete(gamma, endGamma, b)

    # use Scipy.romberg to integrate over spectrum_data_complete between given min and max energy
    prob_sum = romberg(function, minGamma, maxGamma, divmax=20)

    if normalization == 0:
        print("WARNING: normalization value not given, may result in slow performance!")
        print("Calculating normalization...")
        normalization = spectrum_prob_norm(b)

    return prob_sum/normalization


def spectrum_b_range(min_energy, max_energy):
    """
    Calculates spectrum_data_complete function values with
    bFinput and bWMinput set to zero for fitting purposes.
    """
    
    me = 5.10998950e5 #Electron rest mass, in eV.
    endGamma = 3.508e6/me + 1
    minGamma = min_energy/me + 1
    maxGamma = max_energy/me + 1
    maxRange=(3.508e6/me + 1)-(0.0024/me + 1)

    if (maxGamma < minGamma):
        print("Maximum Lorentz gamma must be greater than minimum")
        return 0

    if (max_energy > 3.508e6 or min_energy < 0.0024):
        print("Energies out of range")
        print("Choose energies between 0.0024 - 3.508e6 eV")
        return 0

    def function(gamma):
        return spectrum_data_complete(gamma, endGamma, 0, 0)/gamma

    # use Scipy.romberg to integrate over spectrum_data_complete between given min and max energy
    prob_sum = romberg(function, minGamma, maxGamma, divmax=20)

    return prob_sum

def spec_data_b_derivative(gamma, maxGamma):
    """
    Derivative of spectrum_data_complete with respect to bFinput
    at the specified value of gamma with the specified endpoint gamma
    """
    
    derivative = fc.coulomb_function(3.,gamma,6.) * (maxGamma-gamma)**2 * pow(gamma*gamma-1.0,0.5)
    return derivative
