import csv
import math
import os
import pathlib
import time

import numpy as np
from scipy.interpolate import interp2d
from scipy.misc import derivative


class Field_profile:
    """
    Generates object that calculates field strength of set of coaxial solenoidal coils.
    """

    def __init__(self, list_coils, main_field=0, trap_current=1, interp=True):

        self.trap_current = trap_current
        self._list_coils = list_coils
        self._main_field = main_field

        self._coil_names = []
        self._num_coils = 0

        self._z_offset = 0

        if interp:
            self.initialize_field_strength_interp()

        for coil in list_coils:

            curr_name = coil.name
            if curr_name == None:
                self._num_coils = int(self._num_coils + 1)
                coil.name = "Coil {}".format(self._num_coils)
                self._coil_names.append(curr_name)

            else:
                self._num_coils = int(self._num_coils + 1)
                if not curr_name in self._coil_names:
                    self._coil_names.append(curr_name)
                else:
                    print()
                    print('WARNING: Coil name "{}" already in use'.format(curr_name))
                    print('Using default name "Coil {}"...'.format(self._num_coils))
                    print()
                    coil.name = "Coil {}".format(self._num_coils)
                    self._coil_names.append(curr_name)

    def get_coil(self, name):
        """
        Returns coil with given name.
        """

        for coil in self._list_coils:
            coil_found = False
            if coil.name == name:
                coil_found = True
                return coil
        if coil_found == False:
            print("Coil '{}' does not exist".format(str(name)))
            return None

    def get_coil_list(self):
        """
        Returns list of coils.
        """

        return self._list_coils

    def set_coil_list(self, list_coils):
        """
        Set list of coils.
        """
        self._list_coils = list_coils

    def get_coil_names(self):
        """
        Get names of coils.
        """

        return self._coil_names

    def field_values(self, position, field_coordinates="cylindrical"):
        """
        Calculates the magnetic field T in coordinates field_coordinates at position
        = (radius,phi,z) in cylindrical coordinates with phi in radians. Can choose
        field_coordinates 'Cartesian' or 'cylindrical'.
        """

        Br = 0
        Bz = 0
        for coil in self._list_coils:
            Br = Br + coil.field_values(position, field_coordinates="cylindrical")[0]
            Bz = Bz + coil.field_values(position, field_coordinates="cylindrical")[2]

        Bz = Bz + self._main_field

        if field_coordinates == "cylindrical":

            return (Br, 0, Bz)

        elif field_coordinates == "Cartesian":

            Bx = Br * math.cos(position[1])
            By = Br * math.sin(position[1])

            return (Bx, By, Bz)

        else:
            print("ERROR: {} not a valid coordinate system".format(return_coordinates))

    def field_grad(
        self, position, deriv_order=1, dx=1e-6, grad_coordinates="Cartesian"
    ):

        r = position

        # defining magnetic field component functions
        def Br(rad, ph, z):

            curr_pos = (rad, ph, z)
            return self.field_values(curr_pos, field_coordinates="cylindrical")[0]

        def Bphi(rad, ph, z):

            curr_pos = (rad, ph, z)
            return self.field_values(curr_pos, field_coordinates="cylindrical")[1]

        def Bz(rad, ph, z):

            curr_pos = (rad, ph, z)
            return self.field_values(curr_pos, field_coordinates="cylindrical")[2]

        # calculating gradients
        Bfields = [Br, Bphi, Bz]

        grads = []
        for Bfield in Bfields:

            fr = lambda rad: Bfield(rad, r[1], r[2])
            grad_r = derivative(
                fr, r[0], dx=dx, n=deriv_order, order=(deriv_order * 2 + 1)
            )

            fphi = lambda ph: Bfield(r[0], ph, r[2])
            # grad_phi = derivative(fphi,r[1],dx=dx,n=deriv_order,order=(deriv_order*2 + 1)) / r[0]
            grad_phi = 0

            fz = lambda z: Bfield(r[0], r[1], z)
            grad_z = derivative(
                fz, r[2], dx=dx, n=deriv_order, order=(deriv_order * 2 + 1)
            )

            if grad_coordinates == "cylindrical":
                grads.append((grad_r, grad_phi, grad_z))

            elif grad_coordinates == "Cartesian":

                grad_x = grad_r * math.cos(r[1]) - grad_phi * math.sin(r[1])
                grad_y = grad_r * math.sin(r[1]) + grad_phi * math.cos(r[1])

                grads.append((grad_x, grad_y, grad_z))

        return grads

    def field_strength(self, radius, zpos):
        """
        Calculates Bz magnetic field strength in T at radius, zpos.
        """

        total_field = 0
        for coil in self._list_coils:
            total_field = total_field + coil.field_strength(radius, zpos)

        total_field = total_field + self._main_field
        return total_field

    def initialize_field_strength_interp(self):
        """
        Creates a grid of field strengths based on field_strength() and
        fits it for much faster calls to this function. It also writes
        the calculated grid to a csv for faster excecution.
        """
        waveguide_radius = 0.578e-2  # (m)
        trap_zmax = 5.5e-2  # (m)

        grid_edge_length = 4e-4  # (m), it was found that grid_edge_length = 5e-4 results in 1ppm agreement between field_stength and field_strength_interp

        rho_array = np.arange(0, waveguide_radius, grid_edge_length)
        z_array = np.arange(-trap_zmax, trap_zmax, grid_edge_length)

        dir_path = pathlib.Path(__file__).parent.resolve()
        print(dir_path)
        pkl_path = (
            dir_path
            / "field_profile_pkl_files/main_field_{}_trap_current_{}.csv".format(
                self.main_field, self.trap_current
            )
        )

        try:
            with open(pkl_path, "r") as pkl_file:
                map_array = np.loadtxt(pkl_file)

        except IOError:
            print("Didn't find an existing field map.")
            start = time.process_time()
            map_array = np.zeros((z_array.shape[0], rho_array.shape[0]))

            B = lambda rho, z: self.field_strength(rho, z)

            for i, z in enumerate(z_array):
                for j, rho in enumerate(rho_array):
                    map_array[i][j] = B(rho, z)

            np.savetxt(pkl_path, map_array)
            tot_time = time.process_time() - start
            print("Time to create map_array for new field settings:", tot_time, "\n")
            print("Writing the map_array to csv. \npkl_path: ", pkl_path)

        # Now use the map_array to do the interpolation.
        B_interp2d = interp2d(rho_array, z_array, map_array, kind="cubic")

        # Making it vectorized, meaning float or np array can be inputs.
        def B_interp(rho, z):
            return B_interp2d(rho, z)[0]

        self.field_strength_interpolated_function = np.vectorize(B_interp)

        return None

    def field_strength_interp(self, radius, zpos):
        """
        A faster version of field_strength(), by factor of roughly 1000.
        """

        field_strength = self.field_strength_interpolated_function(radius, zpos)

        return field_strength

    def field_derivative(self, radius, zpos, deriv_order=1, dx=1e-6):
        """
        Calculates dBz / dz derivative on z-axis.
        """

        def f(x):
            return self.field_strength(radius, x)

        return derivative(f, zpos, dx=dx, n=deriv_order, order=(deriv_order * 2 + 1))

    @property
    def main_field(self):
        return self._main_field

    @main_field.setter
    def main_field(self, value):
        self._main_field = value

    @property
    def num_coils(self):
        return self._num_coils

    @property
    def z_offset(self):
        return self._z_offset

    @z_offset.setter
    def z_offset(self, value):

        old_offset = self._z_offset
        for coil in self._list_coils:

            coil.z_center = coil.z_center - old_offset
            coil.z_center = coil.z_center + value
