""" TODO: DOCUMENT

Notes: 
* Get rid of field_strength. Just replace with field_strength_interp.
* Write the field grid to JSON so that it's much faster and more precise.
"""


import os
import json
import math

import pandas as pd
import numpy as np

import he6_cres_spec_sims.spec_tools.spec_calc.spec_calc as sc
import he6_cres_spec_sims.spec_tools.spec_calc.power_calc as pc
from he6_cres_spec_sims.spec_tools.load_default_field_profiles import load_he6_trap

# Math constants.

PI = math.pi
RAD_TO_DEG = 180 / math.pi
P11_PRIME = 1.84118  # First zero of J1 prime (bessel function)

# Physics constants.

ME = 5.10998950e5  # Electron rest mass (eV).
M = 9.1093837015e-31  # Electron rest mass (kg).
Q = 1.602176634e-19  # Electron charge (Coulombs).
C = 299792458  # Speed of light in vacuum, in m/s
J_TO_EV = 6.241509074e18  # Joule-ev conversion


class DotDict(dict):
    """Utility class. dot.notation access to dictionary attributes."""

    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


class Config:
    """TODO: DOCUMENT"""

    def __init__(self, config_filename="config_example"):

        self.load_config_file(config_filename)
        self.load_field_profile()

    def load_config_file(self, config_filename):
        """TODO: DOCUMENT"""

        filepath = "{}/he6_cres_spec_sims/config_files/{}.json".format(
            os.getcwd(), config_filename
        )
        try:
            with open(filepath, "r") as read_file:
                config_dict = json.load(read_file)

                # Take config parameters from config_file.
                self.simulation = DotDict(config_dict["Simulation"])
                self.physics = DotDict(config_dict["Physics"])
                self.hardware = DotDict(config_dict["Hardware"])
                self.kinematics = DotDict(config_dict["Kinematics"])
                self.bandbuilder = DotDict(config_dict["BandBuilder"])
                self.trackbuilder = DotDict(config_dict["TrackBuilder"])

        except Exception as e:
            print("Config file failed to load.")
            raise e

    def load_field_profile(self):
        """TODO: DOCUMENT"""

        try:
            main_field = self.hardware.main_field
            trap_strength = self.hardware.trap_strength

            self.trap_profile = load_he6_trap(main_field, trap_strength)
            self.field_strength = self.trap_profile.field_strength_interp

        except Exception as e:
            print("Field profile failed to load.")
            raise e


class Physics:
    """TODO: DOCUMENT"""

    def __init__(self, config):
        self.config = config

    def generate_beta_energy(self):

        if self.config.physics.monoenergetic == True:
            return self.generate_monoenergetic_beta()

        else:
            raise Exception(
                "Exception: generate_monoenergetic_beta is the only configured setting."
            )

    def generate_monoenergetic_beta(self):
        return self.config.physics.energy

    def generate_beta_position_direction(self):

        position, direction = sc.random_beta_generator(self.config.physics)

        return position, direction


class Hardware:
    """TODO:Document"""

    def __init__(self, config):

        self.config = config
        self.physics = Physics(config)

    def construct_trapped_events_df(self):
        print("~~~~~~~~~~~~Hardware Block~~~~~~~~~~~~~~\n")
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
        initial_zpos = beta_position[2]
        initial_pitch_angle = beta_direction[0]
        initial_phi_dir = beta_direction[1]

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

        if rho_center + max_radius > self.config.hardware.decay_cell_radius:
            print("Not Trapped: Collided with guide wall.")
            trap_condition += 1

        if trap_condition == 0:
            print("Trapped!")
            return True
        else:
            return False


class Kinematics:
    """TODO:Document"""

    def __init__(self, config):

        self.config = config
        self.hardware = Hardware(config)

    def scatter(self, trapped_event_df):
        """TODO: DOCUMENT"""
        print("~~~~~~~~~~~~Kinematics Block~~~~~~~~~~~~~~\n")
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

                # Physics happens. This could maybe be wrapped into a different method.

                # Jump Size: Sampled from normal dist.
                mu = self.config.kinematics.jump_size_eV
                sigma = self.config.kinematics.jump_std_eV
                jump_size_eV = np.random.normal(mu, sigma)

                # Delta Pitch Angle: Sampled from normal dist.
                mu, sigma = 0, self.config.kinematics.pitch_angle_costheta_std
                rand_float = np.random.normal(
                    mu, sigma
                )  # Necessary to properly distribute angles on a sphere.
                delta_center_theta = (np.arccos(rand_float) - PI / 2) * RAD_TO_DEG

                # Second, calculate new pitch angle and energy.
                # New Pitch Angle:
                center_theta = center_theta + delta_center_theta

                # New energy:
                energy = energy_stop - jump_size_eV

                # New position and direction. Only center_theta is changing right now.
                beta_position, beta_direction = (
                    [rho_pos, phi_pos, zpos],
                    [center_theta, phi_dir],
                )

                # Third, construct a scattered, meaning potentially not-trapped, segment df
                scattered_segment_df = self.hardware.construct_untrapped_segment_df(
                    beta_position, beta_direction, energy, event_num
                )

                # Fourth, check to see if the scattered beta is trapped.
                is_trapped = self.hardware.trap_condition(scattered_segment_df)

                jump_num += 1

                # If the event is not trapped or the max number of jumps has been reached,
                # we do not want to write the df to the scattered_segments_list.
                if not is_trapped:
                    print("Event no longer trapped.")
                    break
                if jump_num > self.config.kinematics.jump_num_max:
                    print(
                        "Event reached jump_num_max : {}".format(
                            self.config.kinematics.jump_num_max
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
        main_field = self.config.hardware.main_field
        decay_cell_radius = self.config.hardware.decay_cell_radius

        # Calculate all relevant segment parameters. Order matters here.
        axial_freq = sc.axial_freq(
            df["energy"], df["center_theta"], df["rho_center"], trap_profile
        )
        avg_cycl_freq = sc.avg_cycl_freq(
            df["energy"], df["center_theta"], df["rho_center"], trap_profile
        )
        zmax = sc.max_zpos(df["center_theta"], df["rho_center"], trap_profile)
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
            df["energy"], self.config.hardware.main_field, segment_radiated_power
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
        df["freq_stop"] = freq_stop
        df["energy_stop"] = energy_stop
        df["zmax"] = zmax
        df["mod_index"] = mod_index
        df["slope"] = slope
        df["segment_power"] = segment_power

        return df

    def segment_length(self):
        """TODO: DOCUMENT"""
        mu = self.config.kinematics.mean_track_length
        segment_length = np.random.exponential(mu)

        return segment_length


class BandBuilder:
    """TODO:Document"""

    def __init__(self, config):

        self.config = config

    def bandbuilder(self, segments_df):

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

    def trackbuilder(self, bands_df):

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

        tracks_df = tracks_df.drop(
            columns=[
                "initial_rho_pos",
                "initial_zpos",
                "initial_pitch_angle",
                "trapped_initial_pitch_angle",
                "initial_phi_dir",
                "center_theta",
                "initial_field",
                "initial_radius",
                "center_x",
                "center_y",
                "rho_center",
                "max_radius",
                "zmax",
                "mod_index",
                "avg_cycl_freq",
                "axial_freq",
            ]
        )

        return tracks_df
