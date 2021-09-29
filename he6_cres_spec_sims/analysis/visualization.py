""" TODO: DOCUMENT
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from he6_cres_spec_sims import simulation_blocks as sim_blocks

# Plt settings:
plt.rcParams["figure.dpi"] = 220
plt.rcParams.update({"font.size": 10})


def read_saved_df(config, filename):
    """TODO: DOCUMENT"""

    results_dir = "{}/he6_cres_spec_sims/simulation_results/{}".format(
        os.getcwd(), config.simulation.results_dir
    )
    try:
        df = pd.read_csv(
            "{}/{}.csv".format(results_dir, filename),
            index_col=[0],
        )
    except Exception as e:
        print("File filename not found in {dirpath}.".format(dirpath=results_dir))
        raise e

    return df


def tracks(config, power_color=True, time_interval=1e-8):
    """TODO: DOCUMENT"""

    filename = "trackbuilder_tracks_df"

    # Open tracks_df associated with config.
    try:
        tracks_df = read_saved_df(config, filename)
    except Exception as e:
        print("Have you run the simulation associated with this config?")
        raise e

    # Create visualization of tracks
    run_length = config.trackbuilder.run_length
    power_max = tracks_df["band_power"].max()
    power_min = tracks_df["band_power"].min()

    def power_to_color(power):
        color = (power - power_min) / (power_max - power_min)
        return 1 - color

    for index, row in tracks_df.iterrows():

        time_start = row["time_start"]
        time_stop = row["time_stop"]
        freq_start = row["freq_start"]
        freq_stop = row["freq_stop"]
        slope = (freq_stop - freq_start) / (time_stop - time_start)
        power = row["band_power"]

        time = np.arange(time_start, time_stop, time_interval)
        freq = slope * (time - time_start) + freq_start

        if power_color:
            color = power_to_color(power)
        else:
            color = 0

        plt.plot(time_start, freq_start, "go")
        plt.plot(time, freq, color=str(color))

    plt.figure(1)
    plt.ylabel("Freq (Hz)")
    plt.xlabel("Time (s)")
    plt.title("Simulation Track Visualization")
