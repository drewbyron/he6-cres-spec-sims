import csv
import math
import os
import pathlib
import time

import numpy as np
from scipy.interpolate import interp2d
from scipy.misc import derivative
from scipy.optimize import fmin


class TrapFieldProfile:
    def __init__(self, main_field, trap_current):

        # TODO: May want to protect these variables with underscores?
        # TODO: Would be nice to have an attribute be the relative trap depth.
        # TODO: Add in trap radius as an attribute?
        
        self.trap_current = trap_current
        self.main_field = main_field

        self.field_strength = self.initialize_field_strength_interp()
        self.trap_width = self.trap_width_calc()

        # TODO: Actually test to be sure it is a trap.
        self.is_trap = True

        self.relative_depth = (main_field - self.field_strength(0, 0)) / main_field

    def initialize_field_strength_interp(self):
        """Document"""
        # TODO: hmm I guess these need to be hardcoded for the moment.
        waveguide_radius = 0.578e-2  # (m)
        trap_zmax = 5.5e-2  # (m)

        grid_edge_length = 4e-4  # (m), it was found that grid_edge_length = 5e-4 results in 1ppm agreement between field_stength and field_strength_interp

        rho_array = np.arange(0, waveguide_radius, grid_edge_length)
        z_array = np.arange(-trap_zmax, trap_zmax, grid_edge_length)

        dir_path = pathlib.Path(__file__).parents[0]

        pkl_path = (
            dir_path
            / "trap_field_profile_pkl/2021_trap_profile_mainfield_0T_trap_1A.csv"
        )

        try:
            with open(pkl_path, "r") as pkl_file:
                map_array = np.loadtxt(pkl_file)

        except IOError as e:
            print("Do you have a field map here: {} ".format(pkl_path))
            raise e

        # Adjust the field values so they align with the given trap configuration.
        map_array = map_array * self.trap_current + self.main_field
        # Now use the map_array to do the interpolation.
        B_interp2d = interp2d(rho_array, z_array, map_array, kind="cubic")

        # Making it vectorized, meaning float or np array can be inputs.
        def B_interp(rho, z):
            return B_interp2d(rho, z)[0]

        return np.vectorize(B_interp)

    def trap_width_calc(self):
        """
        Calculates the trap width of the object trap_profile.
        """

        field_func = self.field_strength

        def func(z):
            return -1 * field_func(0, z)

        maximum = fmin(func, 0, xtol=1e-12)[0]
        print("Trap width: ({},{})".format(-maximum, maximum))
        print("Maximum Field: {}".format(-1 * func(maximum)))

        trap_width = (-maximum, maximum)
        return trap_width
