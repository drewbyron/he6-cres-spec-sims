import os

import numpy as np
import pandas as pd
import pytest

from he6_cres_spec_sims import simulation_blocks as sim_blocks

# ~~~~~~~~~~~~~~~~~~~~~~~~Fixtures~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

@pytest.fixture
def config():
    test_config_filename = "test/test_config_0"
    config_instance = sim_blocks.Config(config_filename=test_config_filename)
    return config_instance

@pytest.fixture
def blank_config():
	# Note that main_field and trap_strength can't be blank becuase they
	# are needed to configure a valid trap.
    test_config_filename = "test/test_config_blank"
    config_instance = sim_blocks.Config(config_filename=test_config_filename)
    return config_instance

@pytest.fixture
def blank_band_df():
    blank_band_properties = {
        "energy": 0.0,
        "energy_stop": 0.0,
        "initial_rho_pos": 0.0,
        "initial_phi_pos": 0.0,
        "initial_zpos": 0.0,
        "initial_pitch_angle": 0.0,
        "initial_phi_dir": 0.0,
        "center_theta": 0.0,
        "initial_field": 0.0,
        "initial_radius": 0.0,
        "center_x": 0.0,
        "center_y": 0.0,
        "rho_center": 0.0,
        "trapped_initial_pitch_angle": 0.0,
        "max_radius": 0.0,
        "avg_cycl_freq": 0.0,
        "freq_stop": 0.0,
        "zmax": 0.0,
        "axial_freq": 0.0,
        "mod_index": 0.0,
        "segment_power": 0.0,
        "slope": 0.0,
        "segment_length": 0.0,
        "band_power": 0.0,
        "band_num": 0.0,
        "segment_num": 0,
        "event_num": 0,
    }

    blank_band_df = pd.DataFrame(blank_band_properties, index=[0])
    return blank_band_df
