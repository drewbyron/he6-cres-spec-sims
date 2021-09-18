"""
---
Units for all module functions' inputs/outputs: 

Energy   : eV
B-field  : T
Time     : s
Angle    : degrees
Distance : m
Frequency: Hz
Power    : W


---
"""
# TODO: Organize the imports.
import math
import os

import numpy as np
from numpy.random import uniform

# import json
# from pathlib import Path
# from scipy.integrate import romberg
import scipy.integrate as integrate
from scipy.optimize import fmin
from scipy.optimize import fminbound
from scipy.interpolate import interp1d
import scipy.special as ss

from he6_cres_spec_sims.spec_tools.coil_classes.coil_form import Coil_form
from he6_cres_spec_sims.spec_tools.coil_classes.field_profile import Field_profile
from he6_cres_spec_sims.spec_tools.coil_classes.trap_profile import Trap_profile

# Math constants.

PI = math.pi
RAD_TO_DEG = 180 / math.pi
P11_PRIME = 1.84118  # first zero of J1 prime (bessel functions)

# Physics constants.

ME = 5.10998950e5  # Electron rest mass (eV).
M = 9.1093837015e-31  # Electron rest mass (kg).
Q = 1.602176634e-19  # Electron charge (Coulombs).
C = 299792458  # Speed of light in vacuum, in m/s
J_TO_EV = 6.241509074e18  # Joule-ev conversion


# Simple special relativity functions.


def gamma(energy):
    gamma = (energy + ME) / ME
    return gamma


def momentum(energy):
    momentum = (((energy + ME) ** 2 - ME ** 2) / C ** 2) ** 0.5 / J_TO_EV
    return momentum


def velocity(energy):
    velocity = momentum(energy) / (gamma(energy) * M)
    return velocity


# CRES functions.


def energy_to_freq(energy, field):

    """Converts kinetic energy to cyclotron frequency."""

    cycl_freq = Q * field / (2 * PI * gamma(energy) * M)

    return cycl_freq


def freq_to_energy(frequency, field):

    """Calculates energy of beta particle in eV given cyclotron 
    frequency in Hz, magnetic field in Tesla, and pitch angle 
    at 90 degrees.
    """

    gamma = Q * field / (2 * PI * frequency * M)
    if np.any(gamma < 1):
        gamma = 1
        max_freq = Q * field / (2 * PI * M)
        warning = (
            "Warning: {:.3e} higher than maximum cyclotron frequency {:.3e}".format(
                frequency, max_freq
            )
        )
        print(warning)
    return gamma * ME - ME


def random_beta_generator(parameter_dict):

    """TODO(byron): Think about if the phi_initial parameter has an
    effect.

    Generates a random beta in the trap with pitch angle between
    min_theta and max_theta , and initial position (rho,0,z) between
    min_rho and max_rho and min_z and max_z.
    """

    min_rho = parameter_dict["min_rho"]
    max_rho = parameter_dict["max_rho"]

    min_z = parameter_dict["min_z"]
    max_z = parameter_dict["max_z"]

    min_theta = parameter_dict["min_theta"] / RAD_TO_DEG
    max_theta = parameter_dict["max_theta"] / RAD_TO_DEG

    rho_initial = np.sqrt(uniform(0, 1) * (max_rho ** 2 - min_rho ** 2))
    phi_initial = 2 * PI * uniform(0, 1) * RAD_TO_DEG
    # phi_initial = 0
    z_initial = uniform(min_z, max_z)

    u_min = (1 - np.cos(min_theta)) / 2
    u_max = (1 - np.cos(max_theta)) / 2

    sphere_theta_initial = np.arccos(1 - 2 * (uniform(u_min, u_max))) * RAD_TO_DEG
    sphere_phi_initial = 2 * PI * uniform(0, 1) * RAD_TO_DEG

    position = [rho_initial, phi_initial, z_initial]
    direction = [sphere_theta_initial, sphere_phi_initial]

    return position, direction


def theta_center(zpos, rho, pitch_angle, trap_profile):

    """Calculates the pitch angle an electron with current z-coordinate
    zpos, rho, and current pitch angle pitch_angle takes at the center
    of given trap.
    """

    if trap_profile.is_trap:

        Bmin = trap_profile.field_strength_interp(rho, 0.0)
        Bcurr = trap_profile.field_strength_interp(rho, zpos)

        theta_center_calc = (
            np.arcsin((np.sqrt(Bmin / Bcurr)) * np.sin(pitch_angle / RAD_TO_DEG))
            * RAD_TO_DEG
        )

        return theta_center_calc

    else:
        print("ERROR: Given trap profile is not a valid trap")
        return False


def cyc_radius(energy, field, pitch_angle):

    """Calculates the instantaneous cyclotron radius of a beta electron
    given the energy, magnetic field, and current pitch angle.
    """

    vel_perp = velocity(energy) * np.sin(pitch_angle / RAD_TO_DEG)

    cyc_radius = (gamma(energy) * M * vel_perp) / (Q * field)

    return cyc_radius


def max_radius(energy, center_pitch_angle, rho, trap_profile):

    """Calculates the maximum cyclotron radius of a beta electron given
    thekinetic energy, trap_profile, and center pitch angle (pitch angle
    at center of trap).
    """

    if trap_profile.is_trap:

        min_field = trap_profile.field_strength_interp(rho, 0)
        max_field = min_field / (np.sin(center_pitch_angle / RAD_TO_DEG)) ** 2

        center_radius = cyc_radius(energy, min_field, center_pitch_angle)
        end_radius = cyc_radius(energy, max_field, pitch_angle=90)

        if np.any(center_radius >= end_radius):
            return center_radius
        else:
            print(
                "Warning: max_radius is occuring at end of trap (theta=90). \
                Something odd may be going on."
            )
            return end_radius

    else:
        print("ERROR: Given trap profile is not a valid trap")
        return False


def min_theta(rho, zpos, trap_profile):

    """Calculates the minimum pitch angle theta at which an electron at
    zpos is trapped given trap_profile.
    """

    if trap_profile.is_trap:

        Bmax = trap_profile.field_strength_interp(rho, trap_profile.trap_width[1])
        Bz = trap_profile.field_strength_interp(rho, zpos)

        theta = np.arcsin((Bz / Bmax) ** 0.5) * RAD_TO_DEG
        return theta

    else:
        print("ERROR: Given trap profile is not a valid trap")
        return False


def max_zpos_not_vectorized(center_pitch_angle, rho, trap_profile, debug=False):

    """Calculates the maximum axial length from center of trap as a
    function of center_pitch_angle and rho. Not vectorized. See below
    for vectorized version.
    """

    if trap_profile.is_trap:

        if center_pitch_angle < min_theta(rho, 0, trap_profile):
            print("WARNING: Electron not trapped")
            return False

        else:

            min_field = trap_profile.field_strength_interp(rho, 0)
            max_field = trap_profile.field_strength_interp(
                rho, trap_profile.trap_width[1]
            )

            max_reached_field = min_field / pow(
                math.sin(center_pitch_angle * math.pi / 180), 2
            )

            # initial guess based on treating magnetic well as v-shaped
            slope = (max_field - min_field) / (trap_profile.trap_width[1])
            initial_z = (max_reached_field - min_field) / slope

            if initial_z == trap_profile.trap_width[1]:
                return initial_z

            def func(z):
                curr_field = trap_profile.field_strength_interp(rho, z)
                return abs(curr_field - max_reached_field)

            max_z = fminbound(func, 0, trap_profile.trap_width[1], xtol=1e-14)
            curr_field = trap_profile.field_strength_interp(rho, max_z)

            if (curr_field > max_reached_field) and debug == True:
                print(
                    "Final field greater than max allowed field by: ",
                    curr_field - max_reached_field,
                )
                print("Bmax reached: ", curr_field)

            if debug == True:
                print("zlength: ", max_z)

            if max_z > trap_profile.trap_width[1]:
                print("Error Rate: ", max_z - trap_profile.trap_width[1])

            return max_z

    else:
        print("ERROR: Given trap profile is not a valid trap")
        return False


def max_zpos(center_pitch_angle, rho, trap_profile):

    """Vectorized version of max_zpos_not_vectorized function."""

    max_zpos_vectorized = np.vectorize(max_zpos_not_vectorized)

    return max_zpos_vectorized(center_pitch_angle, rho, trap_profile)


def mod_index(avg_cycl_freq, zmax):

    """Calculates modulation index from average cyclotron frequency
    (avg_cycl_freq) and maximum axial amplitude (zmax).
    """

    # fixed experiment parameters
    waveguide_radius = 0.578e-2
    kc = P11_PRIME / waveguide_radius

    # calculated parameters
    omega = 2 * PI * avg_cycl_freq
    k_wave = omega / C
    beta = np.sqrt(k_wave ** 2 - kc ** 2)

    phase_vel = omega / beta

    mod_index = omega * zmax / phase_vel

    return mod_index


def df_dt(energy, field, power):

    """Calculates cyclotron frequency rate of change of electron with
    given kinetic energy at field in T radiating energy at rate power.
    """

    energy_Joules = (energy + ME) / J_TO_EV

    slope = (Q * field * C ** 2) / (2 * PI) * (power) / (energy_Joules) ** 2

    return slope


def curr_pitch_angle(rho, zpos, center_pitch_angle, trap_profile):

    """Calculates the current pitch angle of an electron at zpos given
    center pitch angle and main field strength.
    """

    if trap_profile.is_trap:

        min_field = trap_profile.field_strength_interp(rho, 0)
        max_z = max_zpos(center_pitch_angle, rho, trap_profile)
        max_reached_field = trap_profile.field_strength_interp(rho, max_z)

        if np.any(abs(zpos) > max_z):
            print("Electron does not reach given zpos")
            curr_pitch = "FAIL"
        else:
            curr_field = trap_profile.field_strength_interp(rho, zpos)
            curr_pitch = np.arcsin(np.sqrt(curr_field / max_reached_field)) * RAD_TO_DEG

        return curr_pitch

    else:
        print("ERROR: Given trap profile is not a valid trap")
        return False


def axial_freq_not_vect(energy, center_pitch_angle, rho, trap_profile):

    """Caculates the axial frequency of a trapped electron."""

    if trap_profile.is_trap:
        # center_pitch_angle = np.where(center_pitch_angle == 90.0, 89.999, center_pitch_angle)
        if center_pitch_angle == 90.0:
            center_pitch_angle = 89.999

        zmax = max_zpos(center_pitch_angle, rho, trap_profile)
        B = lambda z: trap_profile.field_strength_interp(rho, z)
        Bmax = trap_profile.field_strength_interp(rho, zmax)

        # Secant of theta as function of z. Use conserved mu to derive.
        sec_theta = lambda z: (1 - B(z) / Bmax) ** (-0.5)
        T_a = (
            4
            / velocity(energy)
            * integrate.quad(sec_theta, 0, zmax, epsrel=10 ** -2)[0]
        )

        axial_frequency = 1 / T_a

        return axial_frequency

    else:
        print("ERROR: Given trap profile is not a valid trap")
        return False


def axial_freq(energy, center_pitch_angle, rho, trap_profile):

    """Vectorized version of axial_freq_not_vectorized function."""

    axial_freq_vect = np.vectorize(axial_freq_not_vect)

    return axial_freq_vect(energy, center_pitch_angle, rho, trap_profile)


def avg_cycl_freq_not_vect(energy, center_pitch_angle, rho, trap_profile):

    """Calculates the average cyclotron frquency of an electron given
    kinetic energy, main field, and center pitch angle.
    Returns 0 if electron is not trapped.
    """

    if trap_profile.is_trap:

        Bmin = trap_profile.field_strength_interp(rho, 0)
        min_trapped_angle = min_theta(rho, 0, trap_profile)

        if center_pitch_angle < min_trapped_angle:
            print("Warning: electron not trapped")
            return False

        if center_pitch_angle == 90.0:
            avg_cyc_freq = energy_to_freq(energy, Bmin)

        else:
            zmax = max_zpos(center_pitch_angle, rho, trap_profile)
            B = lambda z: trap_profile.field_strength_interp(rho, z)
            Bmax = trap_profile.field_strength_interp(rho, zmax)
            integrand = lambda z: B(z) * ((1 - B(z) / Bmax) ** (-0.5))

            ax_freq = axial_freq(energy, center_pitch_angle, rho, trap_profile)
            
            # Similar to axial_freq calculation.
            avg_cyc_freq = (
                4
                * Q
                * ax_freq
                / (2 * PI * momentum(energy))
                * integrate.quad(integrand, 0, zmax, epsrel=10 ** -2)[0]
            )

        return avg_cyc_freq

    else:
        print("ERROR: Given trap profile is not a valid trap")
        return False


def avg_cycl_freq(energy, center_pitch_angle, rho, trap_profile):

    """Vectorized version of avg_cyc_freq_not_vectorized function."""

    avg_cycl_freq_vect = np.vectorize(avg_cycl_freq_not_vect)

    return avg_cycl_freq_vect(energy, center_pitch_angle, rho, trap_profile)


def t_not_vect(energy, zpos, center_pitch_angle, rho, trap_profile):

    """DEBUG(byron): This is returning negative times.
    Caculates the time for electron to travel from z = 0 to zpos.
    """

    if trap_profile.is_trap:

        if center_pitch_angle == 90.0:
            center_pitch_angle = 89.999

        zmax = max_zpos(center_pitch_angle, rho, trap_profile)
        B = lambda z: trap_profile.field_strength_interp(rho, z)
        Bmax = trap_profile.field_strength_interp(rho, zmax)

        # Secant of theta as function of z. Use conserved mu to derive.
        sec_theta = lambda z: (1 - B(z) / Bmax) ** (-0.5)
        t = (
            1
            / velocity(energy)
            * integrate.quad(sec_theta, 0, zpos, epsrel=10 ** -2)[0]
        )

        if zpos > zmax:
            print("ERROR: zpos equal to or larger than zmax.")

        return t

    else:
        print("ERROR: Given trap profile is not a valid trap")
        return False


def t(energy, zpos, center_pitch_angle, rho, trap_profile):

    """Vectorized version of t_not_vect function."""

    t_vect = np.vectorize(t_not_vect)

    return t_vect(energy, zpos, center_pitch_angle, rho, trap_profile)


def sideband_calc(avg_cycl_freq, axial_freq, zmax, num_sidebands=7):

    """Calculates relative magnitudes of num_sidebands sidebands from
    average cyclotron frequency (avg_cycl_freq), axial frequency
    (axial_freq), and maximum axial amplitude (zmax).
    """

    # #Physical and mathematical constants
    # m = 9.1093837015e-31 # Electron rest mass, in kg.
    # c = 299792458 #Speed of light in vacuum, in m/s
    # p11prime = 1.84118 #first zero of J1 prime

    # #fixed experiment parameters
    # waveguide_radius = 0.578e-2
    # kc = p11prime / waveguide_radius

    # fixed experiment parameters
    waveguide_radius = 0.578e-2
    kc = P11_PRIME / waveguide_radius

    # calculated parameters
    omega = 2 * PI * avg_cycl_freq
    k_wave = omega / C
    beta = np.sqrt(k_wave ** 2 - kc ** 2)

    phase_vel = omega / beta

    mod_index = omega * zmax / phase_vel

    # Calculate K factor
    K = 2 * PI * avg_cycl_freq * zmax / phase_vel

    # Calculate list of (frequency, amplitude) of sidebands
    sidebands = []

    for k in range(-num_sidebands, num_sidebands + 1):

        # if axial_freq == "Indexed":
        #     freq = k
        # else:
        #     freq = avg_cycl_freq + k * axial_freq

        freq = avg_cycl_freq + k * axial_freq
        magnitude = abs(ss.jv(k, K))

        pair = (freq, magnitude)
        sidebands.append(pair)

    return sidebands, mod_index
