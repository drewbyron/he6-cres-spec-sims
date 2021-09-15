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
def physics(config):
    # test_config_filename = "test/test_config_0"
    # config = sim_blocks.Config(config_filename=test_config_filename)
    physics_instance = sim_blocks.Physics(config)
    return physics_instance


@pytest.fixture
def hardware(config):
    # test_config_filename = "test/test_config_0"
    # config = sim_blocks.Config(config_filename=test_config_filename)
    hardware_instance = sim_blocks.Hardware(config)
    return hardware_instance
