from spec_tools.spec_calc import spec_calc as sc
from spec_tools.spec_gen import spec_gen as sg
from spec_tools.spec_classes.he6_energy_spec import He6_energy_spec
import numpy as np
import math
from numpy.random import uniform


class He6_random_beta:
    """
    An object that generates a random He6 beta given Fierz coefficient b with maximum energy error of max_error_diff
    """

    def __init__(self,b=0, max_error_diff = 1000):
    
        self._b = b
        self._max_error_diff = max_error_diff
        self._normalization = sg.spectrum_prob_norm(b)
        
        self._he6_spec = He6_energy_spec(step_size = max_error_diff,b=b)
        
        xvals, yvals = self._he6_spec.energy_vals()
        
        probs = []
        for k in range(len(yvals)):
            curr_probs = sum(yvals[:k])
            probs.append(curr_probs)
            
        self._probs = probs
        self._energies = xvals
        
    def generate_random_beta(self):
    
        print("Calculating random energy...")
        
        prob_value = uniform(0,1)
        curr_probs = self._probs[:]
        curr_energies = self._energies[:]
        
        cut_index = int(round(len(curr_probs)/2,0))

        while len(curr_energies) > 2:
        
            if len(curr_energies) == 3:
                if prob_value >= curr_probs[1]:
                    curr_probs = curr_probs[1:]
                    curr_energies = curr_energies[1:]
                else:
                    curr_probs = curr_probs[:2]
                    curr_energies = curr_energies[:2]
            elif prob_value > curr_probs[cut_index]:
                curr_probs = curr_probs[cut_index:]
                curr_energies = curr_energies[cut_index:]
            else:
                curr_probs = curr_probs[:cut_index + 1]
                curr_energies = curr_energies[:cut_index + 1]
                
            if len(curr_energies) == 2:
                if prob_value < curr_probs[0] or prob_value > curr_probs[1]:
                    print("WARNING: Probability selection failed!")
                    input("PAUSE")
                    
                    print("desired prob:",prob_value)
                    print("min prob:",curr_probs[0],"/ max prob:",curr_probs[-1])
                    print("length:",len(curr_probs))
                
            cut_index = int(round(len(curr_probs)/2,0))
            
        energy = uniform(curr_energies[0],curr_energies[1])
        parameter_dict = {
            "min_rho" : 0,
            "max_rho" : 5.78e-3,
            "min_z" : -5e-3,
            "max_z" : 5e-3,
            "min_theta" : 0,
            "max_theta" : 180
        }
        
        position,direction = sc.random_beta_generator(parameter_dict)
        
        #randomize cylindrical phi position coordinate
        position[1] = 2*math.pi * uniform(0,1) * 180/math.pi
        
        return energy,position,direction
        
    @property
    def b(self):
        return self._b
        
    @property
    def max_error_diff(self):
        return self._max_error_diff
        
        
