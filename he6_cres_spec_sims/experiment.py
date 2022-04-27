from dataclasses import dataclass
from natsort import natsorted
import numpy as np
import pandas as pd
import pathlib
from shutil import copyfile
from shutil import rmtree
import typing
from typing import List
import yaml

from . import simulation as sim
from .simulation_blocks import Config
# from .spec_tools import beta_source as source
import he6_cres_spec_sims.spec_tools.beta_source.beta_source as source
import he6_cres_spec_sims.spec_tools.spec_calc.spec_calc as sc

# Utility function:
def get_experiment_dir(experiment_params: dict) -> pathlib.Path:

    base_config_path = pathlib.Path(experiment_params["base_config_path"])
    experiment_name = experiment_params["experiment_name"]
    parent_dir = base_config_path.parents[0]
    experiment_dir = parent_dir / experiment_name

    return experiment_dir


# Utility function:
def get_config_paths(experiment_params: dict) -> List[pathlib.Path]:

    experiment_dir = get_experiment_dir(experiment_params)

    suffix = "T.yaml"
    config_paths = [
        x
        for x in experiment_dir.glob("**/*{}".format(suffix))
        if (x.is_file() and ("ipynb" not in str(x)))
    ]

    if len(config_paths) == 0:
        raise ValueError(
            "No config files found in experiment_dir: {}. \n\
            Have you run the experiment?\
            (run: exp = exp.Experiment(experiment_params))".format(
                experiment_dir
            )
        )

    # The natsort ensures that you get a good ordering in terms of consecutive integers
    # within the config file name.
    config_paths = natsorted(config_paths, key=str)

    return config_paths


class Experiment:
    def __init__(self, experiment_params: dict) -> None:

        self.experiment_params = experiment_params

        self.run_experiment(self.experiment_params)

    def run_experiment(self, experiment_params: dict) -> None:

        # Make all the config files specified by the experiment dictionary.
        self.create_configs_for_experiment(experiment_params)

        # Make the experiment config file.
        self.create_experiment_config_file(experiment_params)

        # Collect all the config path names.
        self.config_paths = get_config_paths(experiment_params)

        # Run all of those simulations.
        self.run_sims(self.config_paths)

        return None

    def create_configs_for_experiment(self, experiment_params: dict) -> None:

        base_config_path = pathlib.Path(experiment_params["base_config_path"])
        experiment_dir = get_experiment_dir(experiment_params)

        # Make the experiments_dir if it doesn't exist. May want to delete contents here?
        if not experiment_dir.is_dir():
            experiment_dir.mkdir()
            print("Created directory: {} ".format(experiment_dir))

        else:
            print("Directory already exists: {} ".format(experiment_dir))
            print(
                "CAREFUL: Continuing will delete the contents of the above directory.\n"
            )
            input("Press Enter to continue...")
            rmtree(experiment_dir)
            experiment_dir.mkdir()

        # Grab data from experiment_params dict.
        isotope = experiment_params["isotope"]
        events_to_simulate = experiment_params["events_to_simulate"]
        betas_to_simulate = experiment_params["betas_to_simulate"]
        seeds = experiment_params["rand_seeds"]
        fields = experiment_params["fields_T"]
        traps = experiment_params["traps_A"]

        for i, (seed, field, trap) in enumerate(zip(seeds, fields, traps)):

            # Round the field because there are often small rounding errors.
            field = np.around(field, 6)
            config_path = experiment_dir / "{}_field_{}T{}".format(
                i, field, base_config_path.suffix
            )

            copyfile(base_config_path, config_path)

            # Open the config file and grab the contents.
            with open(config_path, "r") as f:
                config_dict = yaml.load(f, Loader=yaml.FullLoader)

            # Make the appropriate altercations to the config_dict
            config_dict["Settings"]["rand_seed"] = int(seed)
            config_dict["Physics"]["events_to_simulate"] = int(events_to_simulate)
            config_dict["Physics"]["betas_to_simulate"] = int(betas_to_simulate)
            config_dict["Physics"]["energy_spectrum"]["beta_source"] = str(isotope)
            config_dict["EventBuilder"]["main_field"] = float(field)
            config_dict["EventBuilder"]["trap_current"] = float(trap)

            with open(config_path, "w") as f:
                yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)

        return None

    def create_experiment_config_file(self, experiment_params: dict) -> None:

        experiment_dir = get_experiment_dir(experiment_params)

        # Define path for experiment config.
        experiment_config = experiment_dir / (
            experiment_params["experiment_name"] + "_exp.yaml"
        )

        with open(experiment_config, "w") as f:
            yaml.dump(experiment_params, f, default_flow_style=False, sort_keys=False)

        return None

    def run_sims(self, config_paths: List[pathlib.Path]) -> None:

        for i, config_path in enumerate(config_paths):
            print("+++++++++++++++++++++++++++++++++++++++++++++++++\n\n")
            print("Running simulation {} / {}\n\n".format(i, len(config_paths)))
            print("+++++++++++++++++++++++++++++++++++++++++++++++++")
            simulation = sim.Simulation(config_path)
            simulation.run_full()

        return None


@dataclass
class ExpResults:

    experiment_params: dict
    base_config: object
    config_paths: List[pathlib.Path]
    sampled_gammas: pd.DataFrame
    experiment_results: pd.DataFrame

    @classmethod
    def load(cls, experiment_params: dict = None, experiment_config_path: str = None, include_sampled_gammas = False):

        # Can either provide the experiment_params dict if you have it or a path to the
        # exp_config.yaml file.

        if experiment_params is None:
            experiment_config_path = pathlib.Path(experiment_config_path)

            # Open the config file and grab the contents.
            with open(experiment_config_path, "r") as f:
                experiment_params = yaml.load(f, Loader=yaml.FullLoader)

        if (experiment_params is None) and (experiment_config_path is None):

            raise ValueError(
                "Need to provide either experiment_params or exp_config_path."
            )

        exp_results_dict = {
            "experiment_params": experiment_params,
            "base_config": Config(experiment_params["base_config_path"]),
            "config_paths": None,
            "sampled_gammas": None,
            "experiment_results": None,
        }

        # Then collect all the config path names.
        config_paths = get_config_paths(experiment_params)
        exp_results_dict["config_paths"] = config_paths

        tracks_list = []
        sampled_gammas = []
        fields = []

        # Figure out how many betas were sampled; depends on the mode (beta_num or event_num).
        # Note that if this is -1 then you get the entire array (event_mode).
        beta_num = experiment_params["betas_to_simulate"]

        for i, config_path in enumerate(config_paths):

            print("+++++++++++++++++++++++++++++++++++++++++++++++++\n\n")
            print("Loading simulation {} / {}\n\n".format(i, len(config_paths)))
            print("+++++++++++++++++++++++++++++++++++++++++++++++++")

            # Get the simulation parameters from the config.
            config = Config(config_path)
            field = config.eventbuilder.main_field
            trap_current = config.eventbuilder.trap_current
            print("\nSet field: {}, Trap current: {}\n".format(field, trap_current))
            results = sim.Results.load(config_path)
            tracks = results.dmtracks
            tracks["simulation_num"] = i
            tracks["field"] = field
            tracks["trap_current"] = trap_current
            tracks_list.append(tracks)

            if include_sampled_gammas: 
                # Get the betas that were sampled during the simulation.
                beta_source_ne19 = source.BetaSource(config)
                sampled_energies = beta_source_ne19.energy_array[: int(beta_num)]
                sampled_gammas.append(sc.gamma(sampled_energies))
                fields.append(field)

        if include_sampled_gammas:
            exp_results_dict["sampled_gammas"] = pd.DataFrame.from_dict(
                dict(zip(fields, sampled_gammas))
            )
        exp_results_dict["experiment_results"] = pd.concat(tracks_list)

        exp_results = cls(
            exp_results_dict["experiment_params"],
            exp_results_dict["base_config"],
            exp_results_dict["config_paths"],
            exp_results_dict["sampled_gammas"],
            exp_results_dict["experiment_results"],
        )

        return exp_results
