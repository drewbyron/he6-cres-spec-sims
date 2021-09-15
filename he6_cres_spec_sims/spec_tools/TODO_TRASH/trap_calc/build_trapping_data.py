import os
import os.path

from spec_tools.spec_calc import spec_calc as sc
import spec_tools.trap_calc.trapping_probability_data_builder as tp
from spec_tools.load_trap_profile import load_trap_profile

def build_trapping_data(config_dict):

    step_size = config_dict["step_size"]
    minFreq = config_dict["minFreq"]
    maxFreq = config_dict["maxFreq"]
    trap_radius = config_dict["trap_radius"]
    field = config_dict["field"]
    trap_profile_dir = config_dict["trap_profile_dir"]
    trap_profile_file = config_dict["trap_profile_file"]
    
    if not trap_profile_dir[-1] == "/":
        trap_profile_dir = trap_profile_dir + "/"
        
    trap_profile = load_trap_profile(trap_profile_dir + trap_profile_file)
    if trap_profile == False:
        print("Trap profile file {} failed to load")
        print("Using default He6-CRES trap profile")
        trap_profile = 0
        
    save_dir = config_dict["save_dir"]
    
    if not save_dir[:-1] == "/":
        save_dir = save_dir + "/"
        
    #Finding repository home directory
    trapping_data_dir  = os.path.dirname(os.path.abspath(__file__))
    trapping_data_dir  = os.path.dirname(trapping_data_dir)

    #path of trapping directory
    trapping_data_dir = os.path.dirname(trapping_data_dir) + "/trapping_data/"

    data_dir = trapping_data_dir + save_dir
    
    if os.path.exists(data_dir) == False:
        print(data_dir + " does not exist, creating...")
        os.mkdir(data_dir)
    
    data_file = "trap_data_{}_{}".format(trap_radius,field)
    
    if os.path.exists(data_dir + data_file + ".json"):
            
        print("File: ",data_file,"exists")
        print("Exiting program...")
        return
    
    #Actual Values
    data_dict = {
        "step_size" : step_size,
        "minFreq" : minFreq,
        "maxFreq" : maxFreq,
        "trap_radii" : [trap_radius],
        "fields" : [field],
        "trap_profile" : trap_profile
    }
    
    tp.save_file("",save_dir + data_file)
    print('File "{}" created'.format(data_file + ".json"))
    
    print(data_file + ".json")
    print("Trap radius:",data_dict["trap_radii"][0])
    print("Trap field:",data_dict["fields"][0])
    
    print("Building trapping probability data...")
    tp.build_data(data_dict,save_dir + data_file)

