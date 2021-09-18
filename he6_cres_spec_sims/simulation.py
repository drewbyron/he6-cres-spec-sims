""" TODO: DOCUMENT
"""


import os
import json
import math

import pandas as pd
import numpy as np

import he6_cres_spec_sims.spec_tools.spec_calc.spec_calc as sc
from he6_cres_spec_sims.spec_tools.load_default_field_profiles import load_he6_trap
from he6_cres_spec_sims import simulation_blocks as sim_blocks


class Simulation:
    """TODO: DOCUMENT"""

    def __init__(self, config):

        self.config = config

    def run(self):

        # Initialize all simulation blocks.
        hardware = sim_blocks.Hardware(self.config)
        kinematics = sim_blocks.Kinematics(self.config)
        bandbuilder = sim_blocks.BandBuilder(self.config)
        trackbuilder = sim_blocks.TrackBuilder(self.config)

        # Create a set of trapped events.
        trapped_events_df = hardware.construct_trapped_events_df()
        self.save_df(trapped_events_df, "hardware_trapped_events_df")

        # Scatter the trapped events, creating segments.
        segments_df = kinematics.scatter(trapped_events_df)
        self.save_df(segments_df, "kinematics_segments_df")

        # Build out the bands of the segments.
        bands_df = bandbuilder.bandbuilder(segments_df)
        self.save_df(bands_df, "bandbuilder_bands_df")

        # Add event start times and drop columns, creating tracks.
        tracks_df = trackbuilder.trackbuilder(bands_df)
        self.save_df(tracks_df, "trackbuilder_tracks_df")

        return None

    def save_df(self, df, filename):
        """TODO: DOCUMENT"""
        print("\n**Block Output:**")

        results_dir = "{}/he6_cres_spec_sims/simulation_results/{}".format(
            os.getcwd(), self.config.simulation.results_dir
        )
        print(
            "{} written to /he6_cres_spec_sims/simulation_results/{}\n".format(
                filename, self.config.simulation.results_dir
            )
        )
        exists = os.path.isdir(results_dir)
        # print(exists)

        # If folder doesn't exist, then create it.
        if not exists:
            os.makedirs(results_dir)
            print("created folder : ", results_dir)

        try:
            df.to_csv("{}/{}.csv".format(results_dir, filename))
        except Exception as e:
            print("Unable to open df.")
            raise e

        return 0
