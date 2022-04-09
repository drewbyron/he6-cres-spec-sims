""" simulation_blocks

This module contains all of the simulation "blocks" used by the 
Simulation class (see simulation.py). Each block simulates the action of
a concrete part of the pipeline from beta creation to a .spec file being
 written to disk by the ROACH DAQ. The module also contains the Config
class that reads the JSON config file and holds all of the configurable
parameters as well as the field profile. An instance of  the Config
class linked to a specific JSON config file is passed to each simulation
block.


The general approach is that pandas dataframes, each row describing a
single CRES data object (event, segment,  band, or track), are passed between
the blocks, each block adding complexity to the simulation. This general
structure is broken by the last two  classes (Daq and SpecBuilder),
which are responsible for creating the .spec (binary) file output. This
.spec file can then be fed into Katydid just as real data would be.

Classes contained in module: 

    * DotDict
    * Config
    * Physics
    * EventBuilder
    * SegmentBuilder
    * BandBuilder
    * TrackBuilder
    * DMTrackBuilder
    * Daq
    * SpecBuilder

"""

import json
import math
import os
import pathlib
import yaml

import numpy as np
from numpy.random import default_rng
import pandas as pd
import matplotlib.pyplot as plt

from he6_cres_spec_sims.daq.frequency_domain_packet import FDpacket
from he6_cres_spec_sims.spec_tools.trap_field_profile import TrapFieldProfile
from he6_cres_spec_sims.spec_tools.beta_source.beta_source import BetaSource
import he6_cres_spec_sims.spec_tools.spec_calc.spec_calc as sc
import he6_cres_spec_sims.spec_tools.spec_calc.power_calc as pc

# TODO: Make the seed a config parameter, and pass rng(seed) around.

rng = default_rng()

# Math constants.

PI = math.pi
RAD_TO_DEG = 180 / math.pi
P11_PRIME = 1.84118  # First zero of J1 prime (bessel function)

# Physics constants.

ME = 5.10998950e5  # Electron rest mass (eV).
M = 9.1093837015e-31  # Electron rest mass (kg).
Q = 1.602176634e-19  # Electron charge (Coulombs).
C = 299792458  # Speed of light in vacuum (m/s)
J_TO_EV = 6.241509074e18  # Joule-ev conversion


class DotDict(dict):
    """Provides dot.notation access to dictionary attributes."""

    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


class Config:
    """A class used to contain the field map and configurable parameters
    associated with a given simulation configuration file (for example:
    config_example.json).

    ...

    Attributes
    ----------
    simulation, physics, eventbuilder, ... : DotDict
        A dictionary containing the configurable parameters associated
        with a given simulation block. The parameters can be accessed
        with dot.notation. For example eventbuilder.main_field would
        return a field value in T.

    trap_profile: Trap_profile
        An instance of a Trap_profile that corresponds to the main_field
        and trap_strength specified in the config file. Many of the
        spec_tool.spec_calc functions take the trap_profile as a
        parameter.

    field_strength: Trap_profile instance method
        Quick access to field strength values. field_strength(rho,z)=
        field magnitude in T at position (rho,z). Note that there is no
        field variation in phi. num_legs : int The number of legs the
        animal has (default 4).

    Methods
    -------
    load_config_file(config_filename)
        Loads the config file.

    load_field_profile()
        Loads the field profile.
    """

    def __init__(self, config_path):
        """
        Parameters
        ----------
        config_filename: str
            The name of the config file contained in the
            he6_cres_spec_sims/config_files directory.
        """
        self.load_config_file(config_path)
        self.load_field_profile()

    def load_config_file(self, config_path):
        """Loads the YAML config file and creates attributes associated
        with all configurable parameters.

        Parameters
        ----------
        config_filename: str
            The name of the config file contained in the
            he6_cres_spec_sims/config_files directory.

        Raises
        ------
        Exception
            If config file isn't found or can't be opened.
        """

        try:
            with open(config_path, "r") as read_file:
                config_dict = yaml.load(read_file, Loader=yaml.FullLoader)

                # Take config parameters from config_file.
                self.physics = DotDict(config_dict["Physics"])
                self.eventbuilder = DotDict(config_dict["EventBuilder"])
                self.segmentbuilder = DotDict(config_dict["SegmentBuilder"])
                self.bandbuilder = DotDict(config_dict["BandBuilder"])
                self.trackbuilder = DotDict(config_dict["TrackBuilder"])
                self.downmixer = DotDict(config_dict["DMTrackBuilder"])
                self.daq = DotDict(config_dict["Daq"])
                self.specbuilder = DotDict(config_dict["SpecBuilder"])

        except Exception as e:
            print("Config file failed to load.")
            raise e

    def load_field_profile(self):
        """Uses the he6 trap geometry (2021), specified in the
        load_he6_trap module, and the main_field and trap strength
        specified in the config file to create an instance of
        Trap_profile.

        Parameters
        ----------
        None

        Raises
        ------
        Exception
            If field profile fails to load.
        """

        try:
            main_field = self.eventbuilder.main_field
            trap_current = self.eventbuilder.trap_current
            # self.trap_profile = load_he6_trap(main_field, trap_current)
            # self.field_strength = self.trap_profile.field_strength_interp

            self.trap_profile = TrapFieldProfile(main_field, trap_current)
            self.field_strength = self.trap_profile.field_strength

        except Exception as e:
            print("Field profile failed to load.")
            raise e


class Physics:
    """TODO: DOCUMENT"""

    def __init__(self, config):

        self.config = config
        self.bs = BetaSource(config)

    def generate_beta_energy(self):

        # Idea: 

        return self.bs.get_energy()
        # if self.config.physics.monoenergetic == True:

        #     return self.generate_monoenergetic_beta()

        # else:

        #     raise NotImplementedError("Only monoenergetic betas have been implemented.")

    # def generate_monoenergetic_beta(self):

    #     return self.config.physics.energy

    def generate_beta_position_direction(self):

        position, direction = sc.random_beta_generator(self.config.physics)

        return position, direction


class EventBuilder:
    """TODO:Document"""

    def __init__(self, config):

        self.config = config
        self.physics = Physics(config)

    def run(self):

        print("~~~~~~~~~~~~EventBuilder Block~~~~~~~~~~~~~~\n")
        print("Constructing a set of trapped events:")
        event_num = 0
        events_to_simulate = self.config.physics.events_to_simulate

        while event_num < events_to_simulate:

            print("\nEvent: {}/{}...\n".format(event_num, events_to_simulate - 1))

            # generate trapped beta
            is_trapped = False

            while not is_trapped:

                (
                    initial_position,
                    initial_direction,
                ) = self.physics.generate_beta_position_direction()

                energy = self.physics.generate_beta_energy()

                single_segment_df = self.construct_untrapped_segment_df(
                    initial_position, initial_direction, energy, event_num
                )

                is_trapped = self.trap_condition(single_segment_df)

            if event_num == 0:

                trapped_event_df = single_segment_df

            else:

                trapped_event_df = trapped_event_df.append(
                    single_segment_df, ignore_index=True
                )

            event_num += 1
        return trapped_event_df

    def construct_untrapped_segment_df(
        self, beta_position, beta_direction, beta_energy, event_num
    ):
        """TODO:Document"""
        # Initial beta position and direction.
        initial_rho_pos = beta_position[0]
        initial_phi_pos = beta_position[1]
        
        initial_pitch_angle = beta_direction[0]
        initial_phi_dir = beta_direction[1]
        initial_zpos = beta_position[2]

        initial_field = self.config.field_strength(initial_rho_pos, initial_zpos)
        initial_radius = sc.cyc_radius(beta_energy, initial_field, initial_pitch_angle)

        # TODO: These center_x and center_y are a bit confusing. May
        # be nice to just build this into the power calc.
        center_x = initial_rho_pos - initial_radius * np.cos(
            (90 - initial_phi_dir) / RAD_TO_DEG
        )
        center_y = initial_radius * np.sin((90 - initial_phi_dir) / RAD_TO_DEG)

        rho_center = np.sqrt(center_x ** 2 + center_y ** 2)

        center_theta = sc.theta_center(
            initial_zpos, rho_center, initial_pitch_angle, self.config.trap_profile
        )

        # Use trapped_initial_pitch_angle to determine if trapped.
        trapped_initial_pitch_angle = sc.min_theta(
            rho_center, initial_zpos, self.config.trap_profile
        )
        max_radius = sc.max_radius(
            beta_energy, center_theta, rho_center, self.config.trap_profile
        )

        segment_properties = {
            "energy": beta_energy,
            "energy_stop": 0.0,
            "initial_rho_pos": initial_rho_pos,
            "initial_phi_pos": initial_phi_pos,
            "initial_zpos": initial_zpos,
            "initial_pitch_angle": initial_pitch_angle,
            "initial_phi_dir": initial_phi_dir,
            "center_theta": center_theta,
            "initial_field": initial_field,
            "initial_radius": initial_radius,
            "center_x": center_x,
            "center_y": center_y,
            "rho_center": rho_center,
            "trapped_initial_pitch_angle": trapped_initial_pitch_angle,
            "max_radius": max_radius,
            "avg_cycl_freq": 0.0,
            "b_avg" : 0.0,
            "freq_stop": 0.0,
            "zmax": 0.0,
            "axial_freq": 0.0,
            "mod_index": 0.0,
            "segment_power": 0.0,
            "slope": 0.0,
            "segment_length": 0.0,
            "band_power": np.NaN,
            "band_num": np.NaN,
            "segment_num": 0,
            "event_num": event_num,
        }

        segment_df = pd.DataFrame(segment_properties, index=[event_num])

        return segment_df

    def trap_condition(self, segment_df):
        """TODO:Document"""
        segment_df = segment_df.reset_index(drop=True)

        if segment_df.shape[0] != 1:
            raise ValueError("trap_condition(): Input segment not a single row.")

        initial_pitch_angle = segment_df["initial_pitch_angle"][0]
        trapped_initial_pitch_angle = segment_df["trapped_initial_pitch_angle"][0]
        rho_center = segment_df["rho_center"][0]
        max_radius = segment_df["max_radius"][0]

        trap_condition = 0

        if initial_pitch_angle < trapped_initial_pitch_angle:
            print("Not Trapped: Pitch angle too small.")
            trap_condition += 1

        if rho_center + max_radius > self.config.eventbuilder.decay_cell_radius:
            print("Not Trapped: Collided with guide wall.")
            trap_condition += 1

        if trap_condition == 0:
            print("Trapped!")
            return True
        else:
            return False


class SegmentBuilder:
    """TODO:Document"""

    def __init__(self, config):

        self.config = config
        self.eventbuilder = EventBuilder(config)

    def run(self, trapped_event_df):
        """TODO: DOCUMENT"""
        print("~~~~~~~~~~~~SegmentBuilder Block~~~~~~~~~~~~~~\n")
        # Empty list to be filled with scattered segments.
        scattered_segments_list = []

        for event_index, event in trapped_event_df.iterrows():

            print("\nScattering Event :", event_index)

            # Assign segment 0 of event with a segment_length.
            event["segment_length"] = self.segment_length()

            # Fill the event with computationally intensive properties.
            event = self.fill_in_properties(event)

            # Extract position and center theta from event.
            center_x, center_y = event["center_x"], event["center_y"]
            rho_pos = event["initial_rho_pos"]
            phi_pos = event["initial_phi_pos"]
            zpos = 0
            center_theta = event["center_theta"]
            phi_dir = event["initial_phi_dir"]

            # Extract necessary parameters from event.
            # TODO(byron): Note that it is slightly incorrect to assume the power doesn't change as time passes.
            energy = event["energy"]
            energy_stop = event["energy_stop"]
            event_num = event["event_num"]

            segment_radiated_power = event["segment_power"] * 2

            # Append segment 0 to scattered_segments_list because segment 0 is trapped by default.
            scattered_segments_list.append(event.values.tolist())

            # Begin with trapped beta (segment 0 of event).
            is_trapped = True
            jump_num = 0

            # The loop breaks when the trap condition is False or the jump_num exceeds self.jump_num_max.
            while True:
                print("Jump: {jump_num}".format(jump_num=jump_num))
                scattered_segment = event.copy()

                # Physics happens. TODO: This could maybe be wrapped into a different method.

                # Jump Size: Sampled from normal dist.
                mu = self.config.segmentbuilder.jump_size_eV
                sigma = self.config.segmentbuilder.jump_std_eV
                jump_size_eV = np.random.normal(mu, sigma)

                # Delta Pitch Angle: Sampled from normal dist.
                mu, sigma = 0, self.config.segmentbuilder.pitch_angle_costheta_std
                rand_float = np.random.normal(
                    mu, sigma
                )  # Necessary to properly distribute angles on a sphere.
                delta_center_theta = (np.arccos(rand_float) - PI / 2) * RAD_TO_DEG

                # Second, calculate new pitch angle and energy.
                # New Pitch Angle:
                center_theta = center_theta + delta_center_theta

                # Solving an issue caused by pitch angles larger than 90.
                if center_theta > 90:
                    center_theta = 180 - center_theta

                # New energy:
                energy = energy_stop - jump_size_eV

                # New position and direction. Only center_theta is changing right now.
                beta_position, beta_direction = (
                    [rho_pos, phi_pos, zpos],
                    [center_theta, phi_dir],
                )

                # Third, construct a scattered, meaning potentially not-trapped, segment df
                scattered_segment_df = self.eventbuilder.construct_untrapped_segment_df(
                    beta_position, beta_direction, energy, event_num
                )

                # Fourth, check to see if the scattered beta is trapped.
                is_trapped = self.eventbuilder.trap_condition(scattered_segment_df)

                jump_num += 1

                # If the event is not trapped or the max number of jumps has been reached,
                # we do not want to write the df to the scattered_segments_list.
                if not is_trapped:
                    print("Event no longer trapped.")
                    break
                if jump_num > self.config.segmentbuilder.jump_num_max:
                    print(
                        "Event reached jump_num_max : {}".format(
                            self.config.segmentbuilder.jump_num_max
                        )
                    )
                    break

                scattered_segment_df["segment_num"] = jump_num
                scattered_segment_df["segment_length"] = self.segment_length()
                scattered_segment_df = self.fill_in_properties(scattered_segment_df)

                scattered_segments_list.append(
                    scattered_segment_df.iloc[0].values.tolist()
                )

                # reset energy_stop, so that the next segment can be scattered based on this energy.
                energy_stop = scattered_segment_df["energy_stop"]

        scattered_df = pd.DataFrame(
            scattered_segments_list, columns=trapped_event_df.columns
        )

        return scattered_df

    def fill_in_properties(self, incomplete_scattered_segments_df):

        """DOCUMENT LATER"""

        df = incomplete_scattered_segments_df.copy()
        trap_profile = self.config.trap_profile
        main_field = self.config.eventbuilder.main_field
        decay_cell_radius = self.config.eventbuilder.decay_cell_radius

        # Calculate all relevant segment parameters. Order matters here.
        axial_freq = sc.axial_freq(
            df["energy"], df["center_theta"], df["rho_center"], trap_profile
        )
        avg_cycl_freq = sc.avg_cycl_freq(
            df["energy"], df["center_theta"], df["rho_center"], trap_profile
        )
        # TODO: Make this more accurate as per discussion with RJ. 
        b_avg = sc.energy_and_freq_to_field(df["energy"], avg_cycl_freq)

        zmax = sc.max_zpos(df["energy"], df["center_theta"], df["rho_center"], trap_profile)
        mod_index = sc.mod_index(avg_cycl_freq, zmax)
        segment_radiated_power = (
            pc.power_calc(
                df["center_x"],
                df["center_y"],
                avg_cycl_freq,
                main_field,
                decay_cell_radius,
            )
            * 2
        )
        slope = sc.df_dt(
            df["energy"], self.config.eventbuilder.main_field, segment_radiated_power
        )

        energy_stop = (
            df["energy"] - segment_radiated_power * df["segment_length"] * J_TO_EV
        )
        freq_stop = sc.avg_cycl_freq(
            energy_stop, df["center_theta"], df["rho_center"], trap_profile
        )
        slope = (freq_stop - avg_cycl_freq) / df["segment_length"]

        segment_power = segment_radiated_power / 2

        df["axial_freq"] = axial_freq
        df["avg_cycl_freq"] = avg_cycl_freq
        df["b_avg"] = b_avg
        df["freq_stop"] = freq_stop
        df["energy_stop"] = energy_stop
        df["zmax"] = zmax
        df["mod_index"] = mod_index
        df["slope"] = slope
        df["segment_power"] = segment_power

        return df

    def segment_length(self):
        """TODO: DOCUMENT"""
        mu = self.config.segmentbuilder.mean_track_length
        segment_length = np.random.exponential(mu)

        return segment_length


class BandBuilder:
    """TODO:Document"""

    def __init__(self, config):

        self.config = config

    def run(self, segments_df):

        print("~~~~~~~~~~~~BandBuilder Block~~~~~~~~~~~~~~\n")
        sideband_num = self.config.bandbuilder.sideband_num
        frac_total_segment_power_cut = (
            self.config.bandbuilder.frac_total_segment_power_cut
        )
        total_band_num = sideband_num * 2 + 1

        band_list = []

        for segment_index, row in segments_df.iterrows():

            sideband_amplitudes = sc.sideband_calc(
                row["avg_cycl_freq"],
                row["axial_freq"],
                row["zmax"],
                num_sidebands=sideband_num,
            )[0]

            for i, band_num in enumerate(range(-sideband_num, sideband_num + 1)):

                if sideband_amplitudes[i][1] < frac_total_segment_power_cut:
                    continue
                else:
                    # copy segment in order to fill in band specific values
                    row_copy = row.copy()

                    # fill in new avg_cycl_freq, band_power, band_num
                    row_copy["avg_cycl_freq"] = sideband_amplitudes[i][0]
                    # Note that the sideband amplitudes need to be squared to give power.
                    row_copy["band_power"] = (
                        sideband_amplitudes[i][1] ** 2 * row.segment_power
                    )
                    row_copy["band_num"] = band_num

                    # append to band_list, as it's better to grow a list than a df
                    band_list.append(row_copy.tolist())

        bands_df = pd.DataFrame(band_list, columns=segments_df.columns)

        return bands_df


class TrackBuilder:
    """TODO:Document"""

    def __init__(self, config):

        self.config = config

    def run(self, bands_df):

        print("~~~~~~~~~~~~TrackBuilder Block~~~~~~~~~~~~~~\n")
        run_length = self.config.trackbuilder.run_length
        events_to_simulate = self.config.physics.events_to_simulate

        # TODO: Event timing is not currently physical.
        # Add time/freq start/stop.
        tracks_df = bands_df.copy()
        tracks_df["time_start"] = np.NaN
        tracks_df["time_stop"] = np.NaN

        tracks_df["freq_start"] = bands_df["avg_cycl_freq"]
        tracks_df["freq_stop"] = (
            bands_df["slope"] * bands_df["segment_length"] + bands_df["avg_cycl_freq"]
        )

        # dealing with timing of the events.
        # for now just put all events in the window... need to think about this.
        trapped_event_start_times = np.random.uniform(0, run_length, events_to_simulate)

        # iterate through the segment zeros and fill in start times.

        for index, row in bands_df[bands_df["segment_num"] == 0.0].iterrows():
            #             print(index)
            event_num = int(tracks_df["event_num"][index])
            #             print(event_num)
            tracks_df["time_start"][index] = trapped_event_start_times[event_num]

        for event in range(0, events_to_simulate):

            # find max segment_num for each event
            segment_num_max = int(
                bands_df[bands_df["event_num"] == event]["segment_num"].max()
            )

            for segment in range(1, segment_num_max + 1):

                fill_condition = (tracks_df["event_num"] == float(event)) & (
                    tracks_df["segment_num"] == segment
                )
                previous_time_condition = (
                    (tracks_df["event_num"] == event)
                    & (tracks_df["segment_num"] == segment - 1)
                    & (tracks_df["band_num"] == 0.0)
                )
                #                 print("previous_time_condition : ", previous_time_condition)
                previous_segment_time_start = tracks_df[previous_time_condition][
                    "time_start"
                ].iloc[0]
                previous_segment_length = tracks_df[previous_time_condition][
                    "segment_length"
                ].iloc[0]

                for index, row in tracks_df[fill_condition].iterrows():
                    tracks_df["time_start"][index] = (
                        previous_segment_time_start + previous_segment_length
                    )

        tracks_df["time_stop"] = tracks_df["time_start"] + tracks_df["segment_length"]

        # tracks_df = tracks_df.drop(
        #     columns=[
        #         "initial_rho_pos",
        #         "initial_zpos",
        #         "initial_pitch_angle",
        #         "trapped_initial_pitch_angle",
        #         "initial_phi_dir",
        #         "center_theta",
        #         "initial_field",
        #         "initial_radius",
        #         "center_x",
        #         "center_y",
        #         "rho_center",
        #         "max_radius",
        #         "zmax",
        #         "mod_index",
        #         "avg_cycl_freq",
        #         "axial_freq",
        #     ]
        # )

        return tracks_df


class DMTrackBuilder:
    """TODO:Document"""

    def __init__(self, config):

        self.config = config

    def run(self, tracks_df):
        """TODO:Document"""

        print("~~~~~~~~~~~~DMTrackBuilder Block~~~~~~~~~~~~~~\n")
        print(
            "DownMixing the cyclotron frequency with a {} GHz signal".format(
                np.around(self.config.downmixer.mixer_freq * 1e-9, 4)
            )
        )
        mixer_freq = self.config.downmixer.mixer_freq

        downmixed_tracks_df = tracks_df.copy()
        downmixed_tracks_df["freq_start"] = (
            downmixed_tracks_df["freq_start"] - mixer_freq
        )
        downmixed_tracks_df["freq_stop"] = downmixed_tracks_df["freq_stop"] - mixer_freq

        return downmixed_tracks_df


class Daq:
    """TODO:Document"""

    def __init__(self, config):

        self.config = config

    def run(self, downmixed_tracks_df):
        """TODO:Document"""
        print("~~~~~~~~~~~~Daq Block~~~~~~~~~~~~~~\n")

        # Allocate power in fW to each bin of the spec_array.
        spec_array = self.allocate_powers_from_df(downmixed_tracks_df)

        return spec_array

    def allocate_powers_from_df(self, downmixed_tracks_df):
        """TODO:Document"""

        daq_freqbw = self.config.daq.daq_freqbw
        freq_bins = self.config.daq.freq_bins
        fft_per_slice = self.config.daq.fft_per_slice
        run_length = self.config.trackbuilder.run_length

        # Start by making the grid in freq and time
        freq_per_bin = daq_freqbw / freq_bins

        time_per_slice = fft_per_slice / freq_per_bin
        time_slices = int(run_length / time_per_slice)

        # Should be almost identical to run_length but a multiple of time_per_slice.
        spec_time_stop = time_slices * time_per_slice

        t_ticks = np.arange(0, spec_time_stop + 1 * time_per_slice, time_per_slice)
        f_ticks = np.arange(0, daq_freqbw + 1 * freq_per_bin, freq_per_bin)

        spec_array = np.zeros((time_slices, freq_bins))
        print("Creating a spec_array with shape: ", spec_array.shape)

        print(
            "\nAllocating power to bins in spec_array for {} tracks.".format(
                downmixed_tracks_df.shape[0]
            )
        )

        # TODO: Would be better to use iterrows.
        for i, dm_track in downmixed_tracks_df.iterrows():
            # for i in range(downmixed_tracks_df.shape[0]):

            time_start, freq_start = (
                dm_track["time_start"],
                dm_track["freq_start"],
            )
            time_stop, freq_stop = (
                dm_track["time_stop"],
                dm_track["freq_stop"],
            )
            t_intercepts = t_ticks[(t_ticks > time_start) & (t_ticks < time_stop)]
            f_intercepts = f_ticks[(f_ticks > freq_start) & (f_ticks < freq_stop)]

            # Find line equation.
            m = (freq_stop - freq_start) / (time_stop - time_start)
            b = freq_start - m * time_start

            # Find list of x values for all grid intersections:
            t_i = (f_intercepts - b) / m
            f_i = t_intercepts * m + b

            t_val_of_grid_intersections = np.concatenate(
                (t_intercepts, t_i, np.asarray([time_start, time_stop]))
            )
            t_val_of_grid_intersections.sort()

            f_val_of_grid_intersections = np.concatenate(
                (f_intercepts, f_i, np.asarray([freq_start, freq_stop]))
            )
            f_val_of_grid_intersections.sort()

            # Now find the midpoints of all the lines:
            t_grid_indx = (
                t_val_of_grid_intersections[:-1] + t_val_of_grid_intersections[1:]
            ) / 2
            f_grid_indx = (
                f_val_of_grid_intersections[:-1] + f_val_of_grid_intersections[1:]
            ) / 2

            # # Now find the indices of the effected cells:
            # t_grid_indx = (t_grid_indx / time_per_slice).astype("int") - 1
            # f_grid_indx = (f_grid_indx / freq_per_bin).astype("int") - 1

            # TODO: RESTORE THIS!
            # Now find the indices of the effected cells. I don't think a minus 1 is necessary here.
            t_grid_indx = (t_grid_indx / time_per_slice).astype("int")
            f_grid_indx = (f_grid_indx / freq_per_bin).astype("int")
            # Find the time portion of the line that traveled through the bin.
            portion_of_tot_band_power = (
                t_val_of_grid_intersections[1:] - t_val_of_grid_intersections[:-1]
            ) / time_per_slice

            # Multiply that portion by the power of that track.
            if self.config.daq.band_power_override:
                #                 print("band_power was overidden with value: ", self.config.daq.band_power_override)
                band_power = self.config.daq.band_power_override
            else:
                band_power = dm_track["band_power"]

            bin_power = portion_of_tot_band_power * band_power

            # Drop indices that are out of the range of the spec file.
            out_of_range_condition = (t_grid_indx < time_slices) & (
                f_grid_indx < freq_bins
            )

            t_grid_indx = t_grid_indx[out_of_range_condition]
            f_grid_indx = f_grid_indx[out_of_range_condition]
            bin_power = bin_power[out_of_range_condition]

            spec_array[t_grid_indx, f_grid_indx] += bin_power

        # Now account for the frequency dependent gain.
        if self.config.daq.gain_override:
            print(
                "\nGain was overidden with value: {} \n".format(
                    self.config.daq.gain_override
                )
            )
            gain = np.full(freq_bins, self.config.daq.gain_override)
        else:
            gain = self.get_measured_gain()

        spec_array *= gain

        # # TODO: DELETE ONCE THIS WORKS:

        # fig, ax = plt.subplots(1, 1, figsize = (14,7))
        # ax.plot(np.nonzero(spec_array)[0]*time_per_slice, np.nonzero(spec_array)[1]*freq_per_bin, "bo", alpha = .5)
        # for index, track in downmixed_tracks_df.iterrows():

        #     TimeCoordinates = (track["time_start"],track["time_stop"])
        #     FreqCoordinates = (track["freq_start"],track["freq_stop"])
        #     ax.plot(TimeCoordinates, FreqCoordinates, 'ro-',markersize=.5,alpha = .5)

        # plt.xlabel("Time (s)")
        # plt.ylabel("Freq (Hz)")
        # plt.title("Visualization of power allocation to spec_array")
        # plt.show()

        # Now account for the freqency dependent noise.
        noise_smooth_1d = self.get_measured_noise()
        noise_array = self.construct_2d_noise_array(time_slices, noise_smooth_1d)

        spec_array += noise_array

        return spec_array

    def get_measured_gain(self):
  
        gain_dir = pathlib.Path(__file__).parents[0] / "daq/gain_noise_measurements/"

        try:
            gain_file_path = gain_dir / "gain.csv"
            gain = np.loadtxt(gain_file_path, delimiter=",")

        except Exception as e:
            print("Unable to open {}/gain.csv".format(gain_dir))
            raise e

        return gain

    def get_measured_noise(self):

        # Now how to open gain and noise:
        noise_dir = pathlib.Path(__file__).parents[0] / "daq/gain_noise_measurements/"

        try:
            noise_file_path = noise_dir / "noise.csv"
            noise_smooth_1d = np.loadtxt(noise_file_path, delimiter=",")

        except Exception as e:
            print("Unable to open {}/noise.csv".format(noise_dir))
            raise e

        return noise_smooth_1d

    def construct_2d_noise_array(self, time_slices, noise_smooth_1d):

        # With chisq dist DOF = 4:
        DOF = 4
        noise_array = np.array(
            [rng.chisquare(DOF, time_slices) * mu / DOF for mu in noise_smooth_1d]
        ).T

        # Recently (10/07/2021) discovered that letting "dtype = uint8"
        # do the rounding is not working well.
        noise_array = np.around(noise_array, 0)

        return noise_array


class SpecBuilder:
    """TODO:Document"""

    def __init__(self, config, config_path):

        self.config_path = config_path
        self.config = config

    def run(self, spec_array):
        """TODO:Document"""
        print("~~~~~~~~~~~~SpecBuilder Block~~~~~~~~~~~~~~\n")

        # Grab the headers for the 4 differnt packets that make up the  (0,1,2,3).
        # TODO: Think of a more elegant way of getting the headers. Opening this file is clunky...
        hdr_list = self.get_headers()

        # Make spec file:
        slices_in_spec = spec_array.shape[0]
        freq_bins = self.config.daq.freq_bins
        freq_bins_per_packet = freq_bins / 4

        specfile_name = self.config.specbuilder.specfile_name

        # First make a results_dir with the same name as the config.
        config_name = self.config_path.stem
        parent_dir = self.config_path.parents[0]
        results_dir = parent_dir / config_name

        exists = results_dir.is_dir()

        # If results_dir doesn't exist, then create it.
        if not exists:
            results_dir.mkdir()
            print("created directory : ", results_dir)

        spec_path = results_dir / "{}.spec".format(specfile_name)

        # Pass "wb" to write a binary file
        with open(spec_path, "wb") as spec_file:

            # Loop over slices in spec file
            for a in range(slices_in_spec * 4):

                packet_num = a % 4
                slice_num = a // 4

                # Write data to out_file
                data = spec_array[slice_num][
                    int(packet_num * freq_bins_per_packet) : int(
                        (packet_num + 1) * freq_bins_per_packet
                    )
                ]
                data = data.astype("uint8")
                # Write appropriate header to spec_file.
                spec_file.write(hdr_list[packet_num])
                # Write data to spec_file.
                data.tofile(spec_file)

        print("\n**Block Output:**")
        print("Successfully wrote a spec file to {} \n".format(spec_path))

        return None

    def get_headers(self):

        header_dir = pathlib.Path(__file__).parents[0]
        header_path = (
            header_dir
            / "daq/example_spec/example_spec_file.spec"
        )

        # open file:
        hdr_list = []
        try:
            with open(header_path, "rb") as in_file:
                for m in range(4):
                    hdr = in_file.read(FDpacket.BYTES_IN_HEADER)
                    hdr_list.append(hdr)
                    data = in_file.read(FDpacket.BYTES_IN_PAYLOAD)

        except Exception as e:
            print(
                "Do you have a roach noise file at {} ?".format(
                    header_path
                )
            )
            raise e

        return hdr_list
