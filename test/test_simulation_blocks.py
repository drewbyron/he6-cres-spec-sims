"""TODO: DOCUMENT"""

import os

import numpy as np
import pandas as pd
import pytest


from he6_cres_spec_sims import simulation_blocks as sim_blocks

# ~~~~~~~~~~~~~~~~~~~~~~~~~Utility~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def create_dataframe(tuple_data):
    """Create pandas df from tuple data with a header."""
    return pd.DataFrame.from_records(tuple_data[1:], columns=tuple_data[0])


# ~~~~~~~~~~~~~~~~~~~~~~~~Tests~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Test Config class.
class TestConfig:
    def test_config_dotdict_0(self, config):
        assert config.physics.events_to_simulate == 2

    def test_config_dotdict_1(self, config):
        assert config.hardware.main_field == 0.689

    def test_config_field_strength_0(self, config):
        center_field = np.around(config.field_strength(0, 0), 6)
        assert center_field == np.around(0.6883110000147292, 6)

    @pytest.mark.skip(reason="test shell")
    @pytest.mark.slow
    def test_field_strength_1(self, config):
        # Compare the interp'd field to the regular field in a way that
        # doesn't take too much time.
        pass


# Test Physics Block.
class TestPhysics:
    def test_Physics(self, config):
        physics = sim_blocks.Physics(config)
        assert physics.config.physics.events_to_simulate == 2


# Test Hardware Block.
class TestHardware:

    # Preparing input data for testing trap_condition.
    data_0 = {
        "initial_pitch_angle": [90],
        "trapped_initial_pitch_angle": [87],
        "rho_center": [0],
        "max_radius": [0],
    }
    data_1 = {
        "initial_pitch_angle": [86],
        "trapped_initial_pitch_angle": [87],
        "rho_center": [0],
        "max_radius": [0],
    }
    data_2 = {
        "initial_pitch_angle": [90],
        "trapped_initial_pitch_angle": [87],
        "rho_center": [0.005],
        "max_radius": [0.001],
    }

    df_0 = pd.DataFrame(data_0)
    df_1 = pd.DataFrame(data_1)
    df_2 = pd.DataFrame(data_2)

    @pytest.mark.parametrize(
        "test_df, expected_trap_condition", [(df_0, True), (df_1, False), (df_2, False)]
    )
    def test_trap_condition(self, blank_config, test_df, expected_trap_condition):
        
        hardware = sim_blocks.Hardware(blank_config)
        trap_condition = hardware.trap_condition(test_df)
        assert trap_condition == expected_trap_condition

    @pytest.mark.skip(reason="test shell")
    def test_construct_untrapped_segment_df(self):
        pass

class TestKinematics:
    

    # @pytest.mark.skip(reason="test shell")
    def test_Kinematics(self, blank_config, blank_band_df):

        # Alter blank_config as needed for test.
        blank_config.kinematics.mean_track_length = 1e-3
        blank_config.kinematics.jump_num_max = 4
        blank_config.hardware.decay_cell_radius = 1e-2

        # Alter blank_band_df as needed for test.
        blank_band_df["energy"] = 20e3
        blank_band_df["center_theta"] = 90.0

        kinematics = sim_blocks.Kinematics(blank_config)
        scattered_df = kinematics.scatter(blank_band_df)
        assert scattered_df.shape[0] == blank_config.kinematics.jump_num_max + 1


class TestBandBuilder:
    

    # @pytest.mark.skip(reason="test shell")
    def test_BandBuilder(self, blank_config, blank_band_df):

        # Alter blank_config as needed for test.
        blank_config.bandbuilder.sideband_num = 10
        blank_config.bandbuilder.frac_total_segment_power_cut = 0.0

        # Alter blank_band_df as needed for test.
        blank_band_df["avg_cycl_freq"] = 1.862263e+10
        blank_band_df["axial_freq"] = 2.787875e+07
        blank_band_df["zmax"] = 0.008227    
        blank_band_df["segment_power"] = 1e-15


        bandbuilder = sim_blocks.BandBuilder(blank_config)
        bandbuilder_df = bandbuilder.bandbuilder(blank_band_df)
        assert bandbuilder_df.shape[0] == blank_config.bandbuilder.sideband_num*2 + 1
        assert np.allclose(bandbuilder_df["band_power"].sum(),blank_band_df["segment_power"])