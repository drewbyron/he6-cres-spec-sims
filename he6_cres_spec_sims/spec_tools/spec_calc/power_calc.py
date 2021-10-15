import numpy as np
import scipy.special as sp
from scipy.integrate import simps
import math

import he6_cres_spec_sims.spec_tools.spec_calc.spec_calc as sc


def rho_phi(center_x, center_y, time, frequency, field):
    """
    Calculates cylindrical rho and phi positional coordinates of electron undergoing 
    cyclotron motion with center of cyclotron orbit at (center_x, center_y) with 
    frequency in Hz and field in Tesla as a function of time in seconds.
    """
    cyc_rad = sc.cyc_radius(sc.freq_to_energy(frequency, field), field, 90)
    omega = 2 * math.pi * frequency

    xcoord = center_x + cyc_rad * math.cos(omega * time)
    ycoord = center_y + cyc_rad * math.sin(omega * time)

    rho = math.sqrt(xcoord ** 2 + ycoord ** 2)
    phi = math.atan2(ycoord, xcoord)

    return (rho, phi)


def fA(center_x, center_y, time, frequency, field, trap_radius=0.578e-2):

    p11prime = 1.84118
    kc = p11prime / trap_radius

    omega = 2 * math.pi * frequency

    rho, phi = rho_phi(center_x, center_y, time, frequency, field)

    fA = math.cos(phi) * math.sin(phi - omega * time) * sp.j1(kc * rho) / (
        kc * rho
    ) - math.sin(phi) * math.cos(phi - omega * time) * sp.jvp(1, kc * rho)

    return fA


def fB(center_x, center_y, time, frequency, field, trap_radius=0.578e-2):

    p11prime = 1.84118
    kc = p11prime / trap_radius

    omega = 2 * math.pi * frequency

    rho, phi = rho_phi(center_x, center_y, time, frequency, field)

    fB = -math.sin(phi) * math.sin(phi - omega * time) * sp.j1(kc * rho) / (
        kc * rho
    ) - math.cos(phi) * math.cos(phi - omega * time) * sp.jvp(1, kc * rho)

    return fB


def fA_fB_integral_sum(center_x, center_y, frequency, field, trap_radius=0.578e-2):

    omega = 2 * math.pi * frequency

    time_values = np.linspace(0, 2 * math.pi / omega, 1000)

    fA_cos_integrand_values = [
        fA(center_x, center_y, time, frequency, field, trap_radius)
        * math.cos(omega * time)
        * omega
        / (2 * math.pi)
        for time in time_values
    ]
    fA_sin_integrand_values = [
        fA(center_x, center_y, time, frequency, field, trap_radius)
        * math.sin(omega * time)
        * omega
        / (2 * math.pi)
        for time in time_values
    ]

    fB_cos_integrand_values = [
        fB(center_x, center_y, time, frequency, field, trap_radius)
        * math.cos(omega * time)
        * omega
        / (2 * math.pi)
        for time in time_values
    ]
    fB_sin_integrand_values = [
        fB(center_x, center_y, time, frequency, field, trap_radius)
        * math.sin(omega * time)
        * omega
        / (2 * math.pi)
        for time in time_values
    ]

    fA_cos_integral = simps(fA_cos_integrand_values, time_values)
    fA_sin_integral = simps(fA_sin_integrand_values, time_values)

    fB_cos_integral = simps(fB_cos_integrand_values, time_values)
    fB_sin_integral = simps(fB_sin_integrand_values, time_values)

    fA_fB_integral_complete = (
        abs(fA_cos_integral) ** 2
        + abs(fA_sin_integral) ** 2
        + abs(fB_cos_integral) ** 2
        + abs(fB_sin_integral) ** 2
    )

    return fA_fB_integral_complete


def power_calc_not_vect(center_x, center_y, frequency, field, trap_radius):

    """Calculates the average cyclotron radiation power (in one direction) in Watts in the
    TE11 mode of an electron undergoing cyclotron motion in the
    cylindrical waveguide around the point (center_x,center_y) with
    frequency in Hz and field in Tesla. 
    """

    q = 1.602176634e-19  # Electron charge, in Coulombs
    c = 299792458  # Speed of light in vacuum, in m/s
    mu0 = 4 * math.pi * 1e-7
    p11prime = 1.84118

    kc = p11prime / trap_radius
    Rcycl = sc.cyc_radius(sc.freq_to_energy(frequency, field), field, 90)
    
    if Rcycl > trap_radius:
        print("Warning: cyclotron radius greater than trap radius")

    # values in power equation
    omega = 2 * math.pi * frequency
    k = omega / c

    if (pow(k, 2) - pow(kc, 2)) >= 0:
        beta = math.sqrt(pow(k, 2) - pow(kc, 2))
    else:
        print("Warning: frequency below TE11 cutoff")
        return 0

    power_constant_numerator = mu0 * pow(omega, 3) * pow(q, 2) * pow(kc * Rcycl, 2)
    power_constant_denominator = (
        math.pi * beta * (pow(p11prime, 2) - 1) * pow(sp.j1(p11prime), 2)
    )

    power = (
        power_constant_numerator
        / power_constant_denominator
        * fA_fB_integral_sum(center_x, center_y, frequency, field, trap_radius)
    )

    return power

def power_calc(center_x, center_y, frequency, field, trap_radius):

    """Vectorized version of power_calc_not_vect function."""

    power_calc_vect = np.vectorize(power_calc_not_vect)

    return power_calc_vect(center_x, center_y, frequency, field, trap_radius)
