import numpy as np

from copy import deepcopy

import json
import os
import os.path

from spec_tools.spec_calc import spec_calc as sc
from spec_tools.coil_classes.trap_profile import Trap_profile
from spec_tools.load_default_field_profiles import load_he6_trap

#trapping probability prob_data elements format:
#element data in prob_data -> [trap_radius,float_field,list_xyvals]
#list_xvals -> [[xval1,yval1],[xval2,yval2],...]


def build_data(build_dict,filename):
    """
    Given dictionary build_dict with appropriate parameters builds trapping probability data dictionary and saves it as JSON file <filename>.
    """

    data_dict = deepcopy(build_dict)
    save_dict = construct_data(data_dict)
    saved = save_file(save_dict,filename)
    
    if saved == True:
        print('Trapping probability data created and saved in "{}.json"'.format(filename))
    else:
        print("Error: data not saved")

def load_data(filename):
    """
    Loads trapping probability data dictionary from JSON file <filename>.
    """

    #Finding repository home directory
    data_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.dirname(data_dir)
    
    #path of trapping data file
    data_dir = os.path.dirname(data_dir) + "/trapping_data"
    data_file = data_dir + "/{}.json".format(filename)

    if not os.path.exists(data_file):
        print('File "{}" does not exist'.format(data_file))
        return False

    #Loading JSON test file
    with open(data_file,"r") as read_file:
        try:
            the_data = json.load(read_file)
        except:
            print('Loaded file "{}.json" does not contain a valid trapping probability data dictionary'.format(filename))
            return False
    try:
        step_size = the_data["step_size"]
        minFreq = the_data["minFreq"]
        maxFreq = the_data["maxFreq"]
        bins = the_data["bins"]
        fields = the_data["fields"]
        trap_radii = the_data["trap_radii"]
        prob_data = the_data["prob_data"]
    except:
        print('Loaded file "{}.json" does not contain a valid trapping probability data dictionary'.format(filename))
        return False
        
    return the_data
    
def add_data(add_dict,filename):
    """
    Adds trapping probability data to an existing trapping probability data dictionary JSON file <filename>.
    """
    
    loaded_data = load_data(filename)
    data_dict = deepcopy(add_dict)
    
    matched_data = compare_dict(data_dict,loaded_data)
    if matched_data == False:
        print("Adding trapping probabiliy data operation aborted...")
        return False
    
    
    prev_data = loaded_data["prob_data"]
    prev_radii = loaded_data["trap_radii"]
    prev_fields = loaded_data["fields"]
    
    
    new_radii = data_dict["trap_radii"]
    new_fields = data_dict["fields"]
    
    duplicate_fields = set(prev_fields) & set(new_fields)
    duplicate_radii = set(prev_radii) & set(new_radii)
    print("duplicate fields:",duplicate_fields)
    print("duplicate radii:",duplicate_radii)
    
    if (set(new_fields) == duplicate_fields) and (set(new_radii) == duplicate_radii):
        print("No new field/radii combinations, data remaining unchanged")
        return
        
    elif set(new_radii) == duplicate_radii:
        temp_new_fields = [field for field in set(new_fields)-duplicate_fields]
        print("new fields:",temp_new_fields)
        
        data_dict["fields"] = temp_new_fields
        added_data = construct_data(data_dict)
        new_data = added_data["prob_data"]
        
    elif set(new_fields) == duplicate_fields:
        temp_new_radii = [radius for radius in set(new_radii) - duplicate_radii]
        print("new radii:",temp_new_radii)
        
        data_dict["trap_radii"] = temp_new_radii
        added_data = construct_data(data_dict)
        new_data = added_data["prob_data"]
        
    else:
        temp_new_radii = [radius for radius in set(new_radii) - duplicate_radii]
        
        data_dict["trap_radii"] = temp_new_radii
        added_data = construct_data(data_dict)
        new_data = added_data["prob_data"]
        
        data_dict["trap_radii"] = prev_radii
        temp_new_fields = [field for field in set(new_fields) - duplicate_fields]
        data_dict["fields"] = temp_new_fields
        added_data = construct_data(data_dict)["prob_data"]
        
        for data in added_data:
            new_data.append(data)
        
    for data in new_data:
        #print("data section:",data)
        #print()
        if not (data[0] in duplicate_radii and data[1] in duplicate_fields):
            prev_data.append(data)
            
    combined_data = prev_data
    combined_data.sort()
    combined_fields =  set(prev_fields) | set(new_fields)
    combined_fields = [field for field in combined_fields]
    combined_fields.sort()
    combined_radii = set(prev_radii) | set(new_radii)
    combined_radii = [radius for radius in combined_radii]
    combined_radii.sort()
    
    data_dict["bins"] = int((data_dict["maxFreq"] - data_dict["minFreq"]) / data_dict["step_size"]) - 1
    data_dict["prob_data"] = combined_data
    data_dict["fields"] = combined_fields
    data_dict["trap_radii"] = combined_radii
    
    
    saved = save_file(data_dict,filename)
    
    if saved == True:
        print('Trapping probability data created and saved in "{}.json"'.format(filename))
    else:
        print("Error: data not saved")
        
def merge_data(data_dict_list,overwrite=False):
    """
    Given list of trapping probability dictionaries, merges
    data into a single dictionary.
    """
    
    step_size = data_dict_list[0]["step_size"]
    minFreq = data_dict_list[0]["minFreq"]
    maxFreq = data_dict_list[0]["maxFreq"]
    
    
    if len(data_dict_list) < 2:
        print("Need 2 or more dictionaries to merge, aborting merge...")
        return False
    
    #testing to make sure metadata matches
    for dict in data_dict_list[1:]:
        if step_size != dict["step_size"]:
            print("Step sizes don't match, aborting merge...")
            return False
        elif minFreq != dict["minFreq"] or maxFreq != dict["maxFreq"]:
            print("Min and max frequencies don't match, aborting merge...")
            return False
            
    if overwrite == False:
            
        #testing to make sure all field lists match
        test_fields = data_dict_list[0]["fields"]
        for dict in data_dict_list[1:]:
            if test_fields != dict["fields"]:
                print("WARNING: Field lists in 'fields' parameters of test spectra data to be merged not equal")
                print("Aborting merging operation...")
                return False
        
        #testing to make sure trapping radii lists don't overlap
        radii_sets = []
        for dict in data_dict_list:
            radii_sets.append(set(dict["trap_radii"]))
          
            isdisjoint = True
            for k in range(len(radii_sets)-1):
                test_set = radii_sets[k]
                compare_sets = radii_sets[(k+1):]
              
                for com_set in compare_sets:
                    if not test_set.isdisjoint(com_set):
                        print((test_set & com_set))
                        isdisjoint = False
                        break
                if isdisjoint == False:
                    break
                      
        if isdisjoint == False:
            print("WARNING: Trap radii lists in 'trap_radii' parameters of test spectra data to be merged not unique")
            print("Aborting merging operation...")
            return False
            
    all_fields = []
    for dict in data_dict_list:
        fields = dict["fields"]
        for field in fields:
            all_fields.append(field)
            
    all_trap_radii = []
    for dict in data_dict_list:
        trap_radii = dict["trap_radii"]
        for radius in trap_radii:
            all_trap_radii.append(radius)
            
    set_trap_radii = set(all_trap_radii)
    all_trap_radii = []
    for radius in set_trap_radii:
        all_trap_radii.append(radius)
            
    set_fields = set(all_fields)
    all_fields = []
    for field in set_fields:
        all_fields.append(field)
       
    all_fields.sort()
    all_trap_radii.sort()
            
    meta_data = {
        "step_size" : step_size,
        "minFreq" : minFreq,
        "maxFreq" : maxFreq,
        "bins" : int((maxFreq - minFreq) / step_size) - 1,
        "fields" : all_fields,
        "trap_radii" : all_trap_radii
    }
    
    all_prob_data = []
    check_repeat = set()
    for dict in data_dict_list:
        prob_data = dict["prob_data"]
        for data in prob_data:
            curr_radius = data[0]
            curr_field = data[1]
            if (curr_radius,curr_field) in check_repeat:
                print("WARNING: Current data set shares same trap radius and field as already added set")
                print("Not adding current data set to merged data...")
            else:
                check_repeat.add((curr_radius,curr_field))
                all_prob_data.append(data)
    
    all_prob_data.sort()
    
    merge_dict = meta_data
    merge_dict["prob_data"] = all_prob_data
    
    print("Trapping probability data successfully merged...")
    
    return merge_dict
            
def construct_data(data_dict):
    """
    Given dictionary data_dict specifying necessary parameters, builds
    trapping probability data as a python dictionary and returns said
    dictionary.
    """


    step_size = data_dict["step_size"]
    minFreq = data_dict["minFreq"]
    maxFreq = data_dict["maxFreq"]
    bins = int((data_dict["maxFreq"] - data_dict["minFreq"]) / data_dict["step_size"]) - 1
    fields = data_dict["fields"]
    trap_radii = data_dict["trap_radii"]
    trap_profile = data_dict["trap_profile"]
    
    if trap_profile == 0:
        print("No trap profile given, loading default He6-CRES trap...")
        trap_profile = load_he6_trap()
    
    fields.sort()
    trap_radii.sort()

    meta_data = {
        "step_size" : data_dict["step_size"],
        "minFreq" : data_dict["minFreq"],
        "maxFreq" : data_dict["maxFreq"],
        "bins" : int((data_dict["maxFreq"] - data_dict["minFreq"]) / data_dict["step_size"]) - 1,
        "fields" : fields,
        "trap_radii" : trap_radii
    }
    

    #Calculate trapping probabilities for each field strength and trap radius
    full_data = []
    for trap_rad in trap_radii:
        for curr_field in fields:
        
            #adjusting trap profile for given field
            trap_profile.main_field = curr_field
        
            #generating frequencies to calculate probabilities
        
            curr_freqs = [i for i in np.linspace(minFreq,
                maxFreq,
                int((maxFreq-minFreq)/step_size))]
      
            xvals = curr_freqs
            yvals = []
        
            for freq in curr_freqs:
                print("Calculating probability for {} Hz...".format(freq))
                yvals.append(sc.trapping_prob(sc.freq_to_energy(freq,
                    curr_field),trap_rad,trap_profile))
        
            zip_iterator = zip(xvals,yvals)
            trap_prob_data = [values for values in zip_iterator]
            full_data.append([trap_rad,curr_field,trap_prob_data])
            
    return_dict = meta_data
    return_dict["prob_data"] = full_data
    
    return return_dict
    
def save_file(save_dict,filename):
    """
    Given dictionary save_dict saves save_dict as JSON file <filename>.
    """

    #Finding repository home directory
    data_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.dirname(data_dir)
    
    #Saving trapping data to JSON file
    data_dir = os.path.dirname(data_dir) + "/trapping_data"
    data_file = data_dir + "/{}.json".format(filename)
    
    with open(data_file,"w") as write_file:
        json.dump(save_dict,write_file)
    
    return True
    
def compare_dict(data_dict1, data_dict2):
    """
    Compares metadata of test spectra dictionaries data_dict1 and data_dict2
    and sends related error message if metadata does not match.
    """
    
    try:
        step_size = data_dict1["step_size"]
        minFreq = data_dict1["minFreq"]
        maxFreq = data_dict1["maxFreq"]
        
        step_size = data_dict2["step_size"]
        minFreq = data_dict2["minFreq"]
        maxFreq = data_dict2["maxFreq"]

    except TypeError:
        print('WARNING: One of the submitted dictionaries is not a valid trapping probability data dictionary')
        return False
    
    if data_dict1["step_size"] != data_dict2["step_size"]:
        print("WARNING: Step sizes do not match!")
        return False
    elif data_dict1["minFreq"] != data_dict2["minFreq"] or data_dict1["maxFreq"] != data_dict2["maxFreq"]:
        print("WARNING: Min and max frequencies do not match!")
        return False
    else:
        return True
    

    
    
