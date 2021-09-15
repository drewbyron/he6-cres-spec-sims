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
    def test_Physics(self, physics):
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
    def test_trap_condition(self, hardware, test_df, expected_trap_condition):
        # print(hardware.cnfg.events_to_simulate)

        trap_condition = hardware.trap_condition(test_df)
        assert trap_condition == expected_trap_condition

    @pytest.mark.skip(reason="test shell")
    def test_construct_untrapped_segment_df(self):
        pass

# Test Kinematics Block.
class TestKinematics:
    @pytest.mark.skip(reason="test shell")
    def test_Kinematics(self):
        pass
