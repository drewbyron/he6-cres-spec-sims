import os
import os.path

import spec_tools.trap_calc.trapping_probability_data_builder as tp
from spec_tools.trap_calc.trap_formula import Trap_formula

def load_trap_formula(trap_data_file = "trap_data/trap_data_0.00578_2.0.json"):
    """
    Loads trap formula using trap data from given trap data file.
    """
    
    print("Checking trapping probability data...")
    print("Trying to load {}...".format(trap_data_file))
    if trap_data_file[-5:] == ".json":
        trap_data_file = trap_data_file[:-5]
        
    loaded_data = tp.load_data(trap_data_file)
    
    if loaded_data == False:
        print("File," (trap_data_file + ".json"),
            "doesn't exist or doesn't contain trapping probability data")
        return False
    else:
        print("Data file: ",(trap_data_file + ".json"),"exists")
                
    print("Performing trapping formula fit...")
    print("Trap radius for trapping formula fit:", loaded_data["trap_radii"][0])
    
    prob_data = loaded_data["prob_data"][0]
                   
    fit_radius = prob_data[0]
    fit_field = prob_data[1]
                    
    fit_xdata = [pair[0] for pair in prob_data[2]]
    fit_ydata = [pair[1] for pair in prob_data[2]]
                    
    fit_formula = Trap_formula(fit_xdata,fit_ydata,fit_field,fit_radius)
                    
    return fit_formula


    
