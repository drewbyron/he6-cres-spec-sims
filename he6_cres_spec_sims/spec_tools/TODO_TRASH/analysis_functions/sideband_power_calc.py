import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
import math

import os
import os.path
import json

from scipy.stats import norm
from scipy.optimize import fmin

from spec_tools.formatting.reformatter import config_tick_reformatter
from spec_tools.formatting.reformatter import sci_reformatter
from spec_tools.formatting.input_checker import promptUser

import spec_tools.spec_calc.spec_calc as sc
from spec_tools.spec_calc.power_calc import power_calc

from spec_tools.load_default_field_profiles import load_he6_trap

#example simulations parameters
"""
simulation_dict = {
    "min_rho" : 0,
    "max_rho" : 5.78e-3,
    "min_z" : 0,
    "max_z" : 4.5e-2,
    "min_theta" : 87.5,
    "max_theta" : 90,
    "frequency" : 20e9,
    "main_field" : 1,
    "trap_strength" : 1e-3,
    "max_simulations" : 10,
    "base_num_sidebands" : 7,
    "sideband_tolerance" : 0.99,
    "simulation_results_dir" : "power_data/power_test_results_test/",
    "simulation_results_file_prefix" : "power_simulation",
    "calculate_axial_frequencies" : True
}

analysis_dict = {
    "num_bins" : 100,
    "main_field" : 1,
    "frequency" : 20e9,
    "simulation_results_dir" : "power_data/power_test_results_test/",
    "simulation_results_file_prefix" : "power_simulation"
}
"""

def plot_sideband_powers(frequency,main_field,center_pitch_angle,
    orbit_center,trap_profile=0):

    fig, ax = plt.subplots(figsize=(9,6))
    
    if trap_profile == 0:
        trap_profile = load_he6_trap()
        
    trap_profile.main_field = main_field

    energy = sc.freq_to_energy(frequency,field=main_field)
    avg_cycl_freq = sc.avg_cyc_freq(energy,
        center_pitch_angle=center_pitch_angle,
        trap_profile=trap_profile)
    zmax = sc.max_zpos(center_pitch_angle,trap_profile=trap_profile)
    axial_freq = sc.axial_freq(energy,
        center_pitch_angle=center_pitch_angle,
        trap_profile=trap_profile)
            
    sidebands = sc.sideband_calc(avg_cycl_freq,axial_freq,zmax,12)[0]
    
    print("Sum of sideband power amplitudes:",sum([pair[1]**2 for pair in sidebands]))
            
    field_sum = []
    for zpos in np.linspace(-zmax,zmax,20):
        field_sum.append(trap_profile.field_strength(orbit_center,zpos))
                
    avg_field = sum(field_sum) / len(field_sum)
    power = power_calc(0,orbit_center,avg_cycl_freq,field=avg_field)
            
    sideband_powers = [[pair[0],pair[1]**2 * power] for pair in sidebands]
           
    xvals = [pair[0]/1e9 for pair in sideband_powers]
    yvals = [pair[1]/1e-15 for pair in sideband_powers]
        
    tick_reformatter_x = config_tick_reformatter(3)
    tick_reformatter_y = config_tick_reformatter(2)
        
    ax.xaxis.set_major_formatter(tick.FuncFormatter(tick_reformatter_x))
    ax.yaxis.set_major_formatter(tick.FuncFormatter(tick_reformatter_y))
        
    ax.set_xlabel(r'Frequency (GHz)')
    ax.set_ylabel(r'Power (fW)')
        
    ax.bar(xvals, yvals, width=0.02)
        
    fig.suptitle(r"Power in Sidebands for {} GHz $\beta$ at Main Field = {} T with Center Pitch Angle {:.2f}$^\circ$".format(round(frequency/1e9,4),main_field,center_pitch_angle) + "\n" +r"and Orbit Center Radius {:.2f}".format(orbit_center))
        #fig.legend(loc="center right")
    fig.subplots_adjust(left=0.16,bottom=0.12)
    fig.show()
    
def run_power_simulations(simulation_dict,plot_generated_beta_path = False,
    plot_sideband_powers=False):

    #set overall parameters
    frequency = simulation_dict["frequency"]
    main_field = simulation_dict["main_field"]
    trap_strength = simulation_dict["trap_strength"]
    
    calculate_axial_frequencies = simulation_dict["calculate_axial_frequencies"]
    
    simulation_results_dir = simulation_dict["simulation_results_dir"]
    simulation_results_file_prefix = simulation_dict["simulation_results_file_prefix"]
    
    #loading trap profile
    trap_profile = load_he6_trap(main_field,trap_strength)
    field_strength = lambda r,z : trap_profile.field_values((r,0,z))[2]
    
    #making sure "/" is at end of directory
    if not simulation_results_dir[-1] == "/":
        simulation_results_dir = simulation_results_dir + "/"
        
    #Finding repository home directory
    home_dir  = os.path.dirname(os.path.abspath(__file__))
    home_dir  = os.path.dirname(home_dir)
    home_dir = os.path.dirname(home_dir)
    
    #Make sure power_data directory exists
    power_data_path = home_dir + "/power_data/"
    
    if os.path.exists(power_data_path) == False:
        print("'power_data/'" + " does not exist, creating...")
        os.mkdir(power_data_path)
        
    #check to make sure simulation results directory exists
    simulation_dir_path = os.path.join(home_dir,simulation_results_dir)
    if not os.path.exists(simulation_dir_path):
    
        prompt = "Directory '{}' does not exist".format(simulation_results_dir[:-1]) + "\nDo you want to create it?(y/n)"
        choice = promptUser(prompt)
        
        if choice == "y":
            os.mkdir(simulation_dir_path)
            print("Directory '{}' created".format(os.path.dirname(simulation_results_dir)))
        else:
            print("Directory '{}' does not exist, aborting...".format(os.path.dirname(simulation_results_dir)))
            return
            
    #find last simulation run
    file_index = 1
    simulation_results_file = "{}_{}.json".format(simulation_results_file_prefix,file_index)
    file_path = os.path.join(os.getcwd(),simulation_results_dir,simulation_results_file)
    while os.path.exists(file_path):
   
        print("File: ",simulation_results_file,"exists")
        file_index = file_index+1
        simulation_results_file = "{}_{}.json".format(simulation_results_file_prefix,file_index)
        file_path = os.path.join(os.getcwd(),simulation_results_dir,simulation_results_file)
        
    num_simulations = file_index
    max_num_simulations = simulation_dict["max_simulations"]
    
    if (num_simulations > max_num_simulations):
        print("Already reached maximum number of simulations, aborting...")
        
    while num_simulations <= max_num_simulations:
        
        print("Running simulation {}/{}...".format(num_simulations,max_num_simulations))
        #generate trapped beta
        is_trapped = False
        while is_trapped == False:
    
            position, direction = sc.random_beta_generator(simulation_dict) # This takes the min/max rho, z, theta as inputs. 
    
            rho_pos = position[0]
            zpos = position[2] # position[1] is the phi position. Could just get rid of for now. 

    
            pitch_angle = direction[0]
            phi_dir = direction[1]
    
            #initial tests
            #rho_pos = 2e-3
            #zpos = 0
            
            #pitch_angle = 89
            #phi_dir = 90
    
            trapped_theta = sc.min_theta(zpos,trap_profile) # math.asin(math.sqrt(Bz / Bmax)) * 180 / math.pi. Trap profile = 0 implies our "normal" trap. 
            
            if (pitch_angle < trapped_theta):
                continue # This will then go back to while statement and generate a new random position and direction. 
                
            center_theta = sc.theta_center(zpos,pitch_angle,trap_profile)
            energy = sc.freq_to_energy(frequency,field=main_field) # Here I can have the energy come from the physics block. 
            curr_field = field_strength(rho_pos,zpos)
            curr_radius = sc.cyc_radius(energy,curr_field,pitch_angle)
            
            
            center_x = rho_pos - curr_radius * math.cos((90-phi_dir)*math.pi/180)
            center_y = curr_radius * math.sin((90-phi_dir)*math.pi/180)
            
            rho_center = math.sqrt(center_x**2 + center_y**2)
            
            max_radius = sc.max_radius(energy,center_theta,trap_profile)
            if (rho_center + max_radius) > simulation_dict["max_rho"]:
                continue
                
            is_trapped = True
        
        print("Generated trapped beta...")
        #print("position:",position)
        #print("direction:",direction)
        #print("center_x:",center_x)
        #print("center_y", center_y)
        
        momentum_x = rho_pos + 0.3 * 5.78e-3 * math.cos(phi_dir*math.pi/180) # What is this? Just momentum direction I assume? 
        momentum_y = 0.3 * 5.78e-3 * math.sin(phi_dir*math.pi/180)
        
        
        if plot_generated_beta_path == True: # I personally feel that the plot generation should be seperated from all of the 
        
            fig, ax = plt.subplots(figsize=(9,6))
        
            circle = plt.Circle((0,0),5.78e-3,fc="none",ec="black",label="trap walls")
            ax.add_patch(circle)
            circle_rad = plt.Circle((center_x,center_y),
                curr_radius,fc="none",ec="red",label="orbital path")
            ax.add_patch(circle_rad)
        
            ax.plot([rho_pos,center_x],[0,center_y],label="orbital radius")
            ax.plot([rho_pos,momentum_x],[0,momentum_y],
                label="transverse velocity direction",color="purple")
            ax.arrow(rho_pos,0,momentum_x-rho_pos,
                momentum_y,width=3e-5,fc="purple",ec="purple",
                length_includes_head="True")
       
            tick_reformatter_x = config_tick_reformatter(3)
            tick_reformatter_y = config_tick_reformatter(2)
            
            ax.xaxis.set_major_formatter(tick.FuncFormatter(tick_reformatter_x))
            ax.yaxis.set_major_formatter(tick.FuncFormatter(tick_reformatter_y))
            
            ax.set_xlabel(r'$x$-coordinate (m)')
            ax.set_ylabel(r'$y$-coordinate (m)')
        
            ax.legend(loc="upper right")
            
            fig_title = r"Orbit of Generated Trapped {} GHz $\beta$".format(round(frequency/1e9,4))
            
            fig_title = fig_title +  " at Main Field = {} T".format(main_field)
            fig_title = fig_title +  r"with Center Pitch Angle {:.2f}$^\circ$ within Trap".format(center_theta)
            fig_title = (fig_title + "\n"
                + r" at Orbit Center ({}, {})".format(sci_reformatter(center_x,3),
                    sci_reformatter(center_y,3)))
                    
            fig.suptitle(fig_title)
            fig.subplots_adjust(left=0.16,bottom=0.12)
            fig.show()
            
            input("Press ENTER to continue...")
            plt.close()
            
        #calculate sideband powers
        print("Calculating power in sidebands...")
        
        energy = sc.freq_to_energy(frequency,field=main_field)
        avg_cycl_freq = sc.avg_cyc_freq(energy,center_theta,trap_profile)
        zmax = sc.max_zpos(center_theta,trap_profile)
        
        if calculate_axial_frequencies == True:
            axial_freq = sc.axial_freq(energy,center_theta,trap_profile)
        else:
            axial_freq = "Indexed"
        
        base_num_sidebands = simulation_dict["base_num_sidebands"]
        sideband_tolerance = simulation_dict["sideband_tolerance"]
        
        #calculate enough sidebands to guarantee fraction of power in sidebands greater than sideband tolerance
        sidebands, mod_index = sc.sideband_calc(avg_cycl_freq,
            axial_freq,
            zmax,
            base_num_sidebands)
        sideband_power_fraction = sum([pair[1]**2 for pair in sidebands])
        curr_num_sidebands = base_num_sidebands
        
        while sideband_power_fraction < sideband_tolerance:
            
            sidebands, mod_index = sc.sideband_calc(avg_cycl_freq,
                axial_freq,
                zmax,
                curr_num_sidebands)
            sideband_power_fraction = sum([pair[1]**2 for pair in sidebands])
            curr_num_sidebands = curr_num_sidebands + 1
            
        print("Sum of sideband power amplitudes:", sideband_power_fraction)
            
        orbit_center_radius = math.sqrt(center_x**2 + center_y**2) # May be intuitive to put this higher in the script? 
        field_sum = []
        for zpos in np.linspace(-zmax,zmax,20):
            field_sum.append(field_strength(orbit_center_radius,zpos))
                       
        avg_field = sum(field_sum) / len(field_sum)
        power = power_calc(center_x,center_y,avg_cycl_freq,field=avg_field)
                   
        sideband_powers = [[pair[0],pair[1]**2 * power] for pair in sidebands]
        
        if plot_sideband_powers == True:
        
            fig, ax = plt.subplots(figsize=(9,6))
            
            if calculate_axial_frequencies == True:
                xvals = [pair[0]/1e9 for pair in sideband_powers]
                bar_width = 0.002
            else:
                xvals = [pair[0] for pair in sideband_powers]
                bar_width = 0.5
                
            yvals = [pair[1]/1e-15 for pair in sideband_powers]
            
            tick_reformatter_x = config_tick_reformatter(3)
            tick_reformatter_y = config_tick_reformatter(2)
            
            ax.xaxis.set_major_formatter(tick.FuncFormatter(tick_reformatter_x))
            ax.yaxis.set_major_formatter(tick.FuncFormatter(tick_reformatter_y))
            
            ax.set_xlabel(r'Frequency (GHz)')
            ax.set_ylabel(r'Power (fW)')
            
            ax.bar(xvals, yvals, width=bar_width)
            
            fig_title = r"Power in Sidebands for {} GHz $\beta$".format(
                round(frequency/1e9,4))
                
            fig_title = fig_title + " at Main Field = {} T".format(main_field)
            fig_title = fig_title + r" with Center Pitch Angle {:.2f}$^\circ$".format(center_theta)
            fig_title = (fig_title + "\n"
                + r"at Orbit Center ({}, {})".format(sci_reformatter(center_x,3),
                sci_reformatter(center_y,3)))
            
            fig.suptitle(fig_title)
            fig.subplots_adjust(left=0.16,bottom=0.12)
            #fig.legend(loc="center right")
            fig.subplots_adjust(left=0.16,bottom=0.12)
            fig.show()
            
            input("Press ENTER to continue...")
            
            plt.close()
        #open data save file and update with new simulation data
        curr_simulation_dict= {
            "rho_pos" : rho_pos,
            "zpos" : zpos,
            "pitch_angle" : pitch_angle,
            "phi_dir" : phi_dir,
            "main_field" : main_field,
            "starting_field" : curr_field,
            "starting_cyc_rad" : curr_radius,
            "center_pitch_angle" : center_theta,
            "center_orbit" : [center_x,center_y],
            "frequency" : frequency,
            "avg_cycl_freq" : avg_cycl_freq,
            "axial_freq" : axial_freq,
            "zmax" : zmax,
            "avg_field" : avg_field,
            "mod_index" : mod_index,
            "total_power" : power,
            "sideband_powers" : sideband_powers,
        }
        
        file_index = 1
        simulation_results_file = "{}_{}.json".format(
            simulation_results_file_prefix,
            file_index)
        file_path = os.path.join(os.getcwd(),
            simulation_results_dir,
            simulation_results_file)
        while os.path.exists(file_path):
        
            #print("File: ",simulation_results_file,"exists")
            file_index = file_index+1
            simulation_results_file = "{}_{}.json".format(simulation_results_file_prefix,file_index)
            file_path = os.path.join(os.getcwd(),simulation_results_dir,simulation_results_file)
        
        if not os.path.exists(file_path):
            with open(file_path,"w") as write_file:
                json.dump(curr_simulation_dict,write_file)
            print("Simulation saved as {}".format(simulation_results_file))
            print()
            
        num_simulations = num_simulations + 1
        
def analyze_power_simulations(analysis_dict,check_sidebands=False):

    mod_indices = []
    total_powers = []
    sideband_powers_array = []
    
    main_field = analysis_dict["main_field"]
    frequency = analysis_dict["frequency"]
    simulation_results_dir = analysis_dict["simulation_results_dir"]
    simulation_results_file_prefix = analysis_dict["simulation_results_file_prefix"]
    
    #making sure "/" is at end of directory
    if not simulation_results_dir[-1] == "/":
        simulation_results_dir = simulation_results_dir + "/"

    #loading power simulations
    file_index = 1
    simulation_results_file = "{}_{}.json".format(simulation_results_file_prefix,file_index)
    file_path = os.path.join(os.getcwd(),simulation_results_dir,simulation_results_file)
    while os.path.exists(file_path):
      
        print("File: ",simulation_results_file,"exists")
        
        with open(file_path,"r") as read_file:
            curr_simulation_dict = json.load(read_file)
            
        if not curr_simulation_dict["main_field"] == main_field:
            print("Current simulated beta's main field not equal to analysis field")
            print("Aborting...")
            return
            
        if not curr_simulation_dict["frequency"] == frequency:
            print("Current simulated beta's frequency not equal to analysis frequency")
            print("Aborting...")
            return
            
        mod_indices.append(curr_simulation_dict["mod_index"])
        sideband_powers_array.append(curr_simulation_dict["sideband_powers"])
        total_powers.append(curr_simulation_dict["total_power"])
        
        #load next simulation file
        file_index = file_index+1
        simulation_results_file = "{}_{}.json".format(simulation_results_file_prefix,file_index)
        file_path = os.path.join(os.getcwd(),
            simulation_results_dir,
            simulation_results_file)
        
    max_sideband_powers = []
    for sideband_powers in sideband_powers_array:
        
        #print("Current sideband powers:",sideband_powers)
        curr_sideband_powers = [pair[1] for pair in sideband_powers]
        max_sideband_power = max(curr_sideband_powers)
        #print("Maximum sideband power:",max_sideband_power)
        max_sideband_powers.append(max_sideband_power)
    
    if check_sidebands == True:
    
        for k in range(len(sideband_powers_array)):
        
            curr_sideband_powers = sideband_powers_array[k]
            #print(curr_sideband_powers)
            sum_sidebands = sum([pair[1] for pair in curr_sideband_powers]) / total_powers[k]
            print("Sum sideband power amplitudes:",sum_sidebands)
            input("PAUSE")
    
    fig, ax = plt.subplots(figsize=(9,6))
        
    xhist = [x / 1e-15 for x in max_sideband_powers]
    xhist2 = [x / 1e-15 for x in total_powers]
    num_bins = analysis_dict["num_bins"]
            
    tick_reformatter_x = config_tick_reformatter(2)
    ax.xaxis.set_major_formatter(tick.FuncFormatter(tick_reformatter_x))
    ax.set_xlabel(r'Total Power (fW)')
    ax.set_ylabel(r'Counts')
    
    n, bins, patches = ax.hist(xhist,num_bins,facecolor="blue", alpha = 0.5, label="Doppler-shifted")
    n2, bins2, patches2 = ax.hist(xhist2,num_bins,facecolor="red", alpha = 0.5,label="Non-Doppler-shifted")
        
    fig.suptitle(r"Maximum Power for {} GHz $\beta$ at Main Field = {} T".format(round(frequency/1e9,4),main_field))
    fig.subplots_adjust(left=0.16,bottom=0.12)
    ax.legend(loc="upper right")
    fig.subplots_adjust(left=0.16,bottom=0.12)
    fig.show()
    
def analyze_power_slopes(analysis_dict,check_sidebands=False):

    mod_indices = []
    total_powers = []
    sideband_powers_array = []
    
    main_field = analysis_dict["main_field"]
    frequency = analysis_dict["frequency"]
    cutoff_ratio = analysis_dict["cutoff_ratio"]
    simulation_results_dir = analysis_dict["simulation_results_dir"]
    simulation_results_file_prefix = analysis_dict["simulation_results_file_prefix"]
    
    #making sure "/" is at end of directory
    if not simulation_results_dir[-1] == "/":
        simulation_results_dir = simulation_results_dir + "/"

    #loading power simulations
    file_index = 1
    simulation_results_file = "{}_{}.json".format(simulation_results_file_prefix,file_index)
    file_path = os.path.join(os.getcwd(),simulation_results_dir,simulation_results_file)
    while os.path.exists(file_path):
      
        print("File: ",simulation_results_file,"exists")
        
        with open(file_path,"r") as read_file:
            curr_simulation_dict = json.load(read_file)
            
        if not curr_simulation_dict["main_field"] == main_field:
            print("Current simulated beta's main field not equal to analysis field")
            print("Aborting...")
            return
            
        if not curr_simulation_dict["frequency"] == frequency:
            print("Current simulated beta's frequency not equal to analysis frequency")
            print("Aborting...")
            return
            
        mod_indices.append(curr_simulation_dict["mod_index"])
        sideband_powers_array.append(curr_simulation_dict["sideband_powers"])
        total_powers.append(curr_simulation_dict["total_power"])
        
        #load next simulation file
        file_index = file_index+1
        simulation_results_file = "{}_{}.json".format(simulation_results_file_prefix,file_index)
        file_path = os.path.join(os.getcwd(),
            simulation_results_dir,
            simulation_results_file)
        
    slope_powers = []
    for k in range(len(sideband_powers_array)):
    
        sideband_powers = sideband_powers_array[k]
        curr_total_power = total_powers[k]
        
        #calculate index of center of sideband powers
        carrier_index = (len(sideband_powers) // 2)
        
        #curr_sideband_labels = [pair[0] for pair in sideband_powers]
        #carrier_index = curr_sideband_labels.index(0)
       
        carrier_power = sideband_powers[carrier_index][1]
        
        curr_energy = sc.freq_to_energy(frequency, main_field)
        curr_slope = sc.df_dt(curr_energy,main_field,curr_total_power)
       
        slope_powers.append([curr_slope,carrier_power])
        
    if check_sidebands == True:
    
        for k in range(len(sideband_powers_array)):
        
            curr_sideband_powers = sideband_powers_array[k]
            #print(curr_sideband_powers)
            sum_sidebands = (sum([pair[1] for pair in curr_sideband_powers]) /
                total_powers[k])
            print("Sum sideband power amplitudes:",sum_sidebands)
            input("PAUSE")
    
    fig, ax = plt.subplots(figsize=(9,6))
        
    max_power = max([pair[1] for pair in slope_powers])
    
    xvals = []
    yvals = []
    for pair in slope_powers:
    
        if pair[1] >= max_power * cutoff_ratio:
            xvals.append(pair[0])
            yvals.append(pair[1])
        
    yvals = [y/1e-15 for y in yvals]
 
    ax.scatter(xvals,yvals,s=4,alpha=0.5)
            
    tick_reformatter_x = config_tick_reformatter(2)
    ax.xaxis.set_major_formatter(tick.FuncFormatter(tick_reformatter_x))
    ax.set_xlabel(r'Total Power (Slope) (Hz/s)')
    ax.set_ylabel(r'Carrier Power (Detected) (fW)')
    
    fig.suptitle(r"Total Power Slopes for {} GHz $\beta$ at Main Field = {} T".format(round(frequency/1e9,4),main_field))
    fig.subplots_adjust(left=0.16,bottom=0.12)
    #ax.legend(loc="upper right")
    fig.subplots_adjust(left=0.16,bottom=0.12)
    fig.show()
    

