import numpy as np
from scipy import special
from mpmath import fp
import matplotlib.pyplot as plt
from scipy import integrate

# Can we emake these more precise?
alpha = 1.0 / 137.0
me = 511 * 1e3
hbarc = 197.3
mprot = 1836.15267389


def p(W):
    return np.sqrt(W**2 - 1)


def beta(W):
    return p(W) / W


# Z is for residual nucleus
def eta(W, Z):
    return alpha * Z * W / p(W)


def Gamma(Z):
    return np.sqrt(1 - (alpha * Z) ** 2)


def R(A):
    return 1.2 * A ** (1.0 / 3.0) * 1e-6 * me / hbarc


def F(W, Z, A):
    return (
        2
        * (1 + Gamma(Z))
        * (2 * p(W) * R(A)) ** (2 * (Gamma(Z) - 1))
        * np.exp(-np.pi * eta(W, Z))
        * np.abs(
            special.gamma(Gamma(Z) + 1j * eta(W, Z)) / special.gamma(2 * Gamma(Z) + 1)
        )
        ** 2
    )


def L0(W, Z, A):
    r = R(A)
    # breaking it up by powers of alphaZ
    az1 = (
        alpha
        * Z
        * (
            28.0 / 15.0 * W * r
            - 8.0 / 15.0 * r / W
            + ((W * r) ** 3) * (1.0 + 9.0 / 2.0 * alpha * Z)
        )
    )
    az2 = 7.0 / 20.0 * (alpha * Z) ** 2
    az3 = -0.5 * W * r * (alpha * Z) ** 3
    az4 = -3.0 / 8.0 * W * r * (alpha * Z) ** 4
    az6 = (alpha * Z) ** 6 * (1.0 / 20.0 + 10 * (W * r) ** 2)
    az8 = -3.0 / 8.0 * (alpha * Z) ** 8
    return 1.0 + az1 + az2 + az3 + az4 + az6 + az8 - 1.0 / 3.0 * (p(W) * r) ** 2


def Delta(W, Z, A):
    return (2 * L0(W, Z, A)) / (1 + Gamma(Z)) - 1.0


polylog = np.vectorize(fp.polylog)


def g(W, Wmax):
    const = 3.0 * np.log(mprot) - 3.0 / 4.0
    deltaW = (
        4
        * (np.arctanh(beta(W)) / beta(W) - 1.0)
        * ((Wmax - W) / (3 * W) - 3 / 2 + np.log(2 * (Wmax - W)))
    )
    polybeta = -4.0 / beta(W) * polylog(2, (2 * beta(W)) / (1 + beta(W)))
    deltaW2 = (
        1.0
        / beta(W)
        * np.arctan(beta(W))
        * (
            2.0 * (1.0 + (beta(W)) ** 2)
            + (Wmax - W) ** 2 / (6 * W**2)
            - 4 * np.arctan(beta(W))
        )
    )
    return const + deltaW + polybeta + deltaW2


def dNdE0(W, Wmax):
    return p(W) * W * (Wmax - W) ** 2


def dNdEnaive(W, Wmax, Z, A):
    return F(W, Z, A) * dNdE0(W, Wmax)


def dNdE(W, Wmax, Z, A):
    def dNdE_unnormed(W, Wmax, Z, A):
        return (
            dNdEnaive(W, Wmax, Z, A)
            * (1.0 + Delta(W, Z, A))
            * (1.0 + alpha / (2.0 * np.pi))
            * g(W, Wmax)
        )

    norm, norm_err = integrate.quad(
        dNdE_unnormed,
        1.000001,
        Wmax,
        args=(Wmax, Z, A),
    )

    return dNdE_unnormed(W, Wmax, Z, A) / norm


class BetaSpectrum:
    """A class used to contain the field map and configurable parameters
    associated with a given simulation configuration file (for example:
    config_example.json).

    ...

    Attributes
    ----------


    Methods
    -------

    """

    def __init__(self, isotope_info: dict ) -> None:
        """
        Parameters
        ----------
        config_filename: str
            The name of the config file contained in the
            he6_cres_spec_sims/config_files directory.
        """
        self.dNdE_norm = None
        self.Wmax = isotope_info["Wmax"]
        self.Z = isotope_info["Z"]
        self.A = isotope_info["A"]

    def dNdE(self, W):
        def dNdE_unnormed(W, Wmax, Z, A):
            return (
                dNdEnaive(W, Wmax, Z, A)
                * (1.0 + Delta(W, Z, A))
                * (1.0 + alpha / (2.0 * np.pi))
                * g(W, Wmax)
            )

        if self.dNdE_norm is None:

            norm, norm_err = integrate.quad(
                dNdE_unnormed,
                1.0,
                self.Wmax,
                args=(self.Wmax, self.Z, self.A),
            )
            self.dNdE_norm = norm

        return dNdE_unnormed(W, self.Wmax, self.Z, self.A) / self.dNdE_norm
