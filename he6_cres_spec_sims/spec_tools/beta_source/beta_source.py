# Local modules.
from .beta_spectrum import BetaSpectrum


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

        # # Reorganize this!! 
        # self.freq_range = (18.0e9, 19.1e9)

        # Include all the sources here.
        self.isotopes = {
            "Ne19": {"Wmax": 2216400.0, "Z": 9, "A": 19},
            "He6": {"Wmax": 0, "Z": 0, "A": 0},
        }

        if (self.source["beta_source"] == "Ne19") or (
            self.source["beta_source"] == "He6"
        ):

            self.beta_spectrum = BetaSpectrum(self.isotopes[self.source["beta_source"]])

    def get_energy(self):

        # Randomly sample (carefully) and then see if it's in the freq domain you want. 
        # Check that this recreates the spectrum slice by slice before moving on! (by changing the field)
        W = 1.05
        return self.beta_spectrum.dNdE(W)
