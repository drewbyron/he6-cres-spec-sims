import numpy as np

# Local modules.
from .beta_spectrum import BetaSpectrum
from ..spec_calc import spec_calc as sc

ME = 5.10998950e5  # Electron rest mass (eV).


class BetaSource:
    """A class used to ...

    ...

    Attributes
    ----------


    Methods
    -------

    """

    def __init__(self, config):
        """
        DOCUMENT
        """
        self.config = config
        self.source = self.config.physics.energy_spectrum
        self.energy_array_len = 10**6
        self.energy_array = None

        # Include all the sources here.
        self.isotopes = {
            "Ne19": {"Wmax": 4.337377690802349, "Z": 9, "A": 19},
            "He6": {"Wmax": 0, "Z": 0, "A": 0},
        }

        if (self.source["beta_source"] == "Ne19") or (
            self.source["beta_source"] == "He6"
        ):
            print("Ne19 or He6 source.")
            self.beta_spectrum = BetaSpectrum(self.isotopes[self.source["beta_source"]])

        self.make_energy_samples()

    def make_energy_samples(self):

        if (self.source["beta_source"] == "Ne19") or (
            self.source["beta_source"] == "He6"
        ):
            if type(self.energy_array) != np.ndarray:

                print("Creating energy_array to pull beta energies from.")
                # Freq acceptance cut:
                main_field = self.config.eventbuilder.main_field

                freq_acceptance_high = self.config.physics.freq_acceptance_high
                freq_acceptance_low = self.config.physics.freq_acceptance_low

                self.energy_acceptance_high = sc.freq_to_energy(
                    freq_acceptance_high, main_field
                )
                self.energy_acceptance_low = sc.freq_to_energy(
                    freq_acceptance_low, main_field
                )

                # Now energy_array has units of eV.
                (
                    self.energy_array,
                    self.fraction_of_spectrum,
                ) = self.beta_spectrum.energy_samples(
                    self.energy_array_len,
                    self.energy_acceptance_low,
                    self.energy_acceptance_high,
                )
                print((
                    self.energy_array,
                    self.fraction_of_spectrum,
                ))

            return None

        else:

            raise NotImplementedError("Kr needs to be reimplimented.")

        # Be sure this doesn't add any weird sampling issues.
        # This will break if beta_num is larger than energy_array_len. Think about.
        print(self.energy_array[beta_num])
        return self.energy_array[beta_num]
