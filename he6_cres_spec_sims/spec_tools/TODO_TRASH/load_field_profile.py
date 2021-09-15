import json
import os
import os.path

from spec_tools.coil_classes.coil_form import Coil_form
from spec_tools.coil_classes.field_profile import Field_profile

from scipy.misc import derivative
from scipy.optimize import fminbound


def load_field_profile(filename):

    #Finding repository home directory
    home_dir  = os.path.dirname(os.path.abspath(__file__))
    home_dir  = os.path.dirname(home_dir)
    
    
    if not filename[0] == "/":
        data_file = home_dir + "/" + filename
    else:
        data_file = home_dir + filename
 
    #Loading JSON test file
    with open(data_file,"r") as read_file:
        try:
            config_dict = json.load(read_file)
            field_coils = config_dict["field_coils"]
            main_field = config_dict["main_field"]
        except:
            print('Loaded file "{}" does not contain a valid field profile config dictionary'.format(filename))
            return
            
    print("Field profile {} loaded".format(filename))

    list_coils= []
    for coil_dict in field_coils:
        
        inner_radius = coil_dict["inner_radius"]
        outer_radius = coil_dict["outer_radius"]
        left_edge = coil_dict["left_edge"]
        right_edge = coil_dict["right_edge"]
        z_center = coil_dict["z_center"]
        num_windings = coil_dict["num_windings"]
        current_per_wire = coil_dict["current_per_wire"]
        name = coil_dict["name"]

        curr_coil = Coil_form(inner_radius=inner_radius,
        outer_radius=outer_radius,
        left_edge=left_edge,
        right_edge=right_edge,
        z_center=z_center,
        num_windings=num_windings,
        current_per_wire=current_per_wire,
        name=name)
        
        list_coils.append(curr_coil)

    field_profile = Field_profile(list_coils,main_field)

    return field_profile
