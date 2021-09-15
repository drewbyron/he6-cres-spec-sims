from he6_cres_spec_sims.spec_tools.coil_classes.coil_form import Coil_form
from he6_cres_spec_sims.spec_tools.coil_classes.field_profile import Field_profile
from he6_cres_spec_sims.spec_tools.coil_classes.trap_profile import Trap_profile
# from spec_tools.spec_calc.spec_calc import generate_triple_coil_field


from scipy.misc import derivative
from scipy.optimize import fminbound


def load_main_magnet(main_field=1):

    current_per_wire = 88 / 7 * main_field
    
    center_coils_radii = [1.*12.7e-2 , 3*12.7e-2]
    inner_coils_radii = [1.25*12.7e-2 , 1.35*12.7e-2]
    mid_coils_radii = [1*12.7e-2 , 3.75*12.7e-2]
    outer_coils_radii = [1*12.7e-2 , 1.15*12.7e-2]
    trim_coils_radii = [1*0.127 , 3*0.127]

    center_coils_edge = 11e-2
    inner_coils_edge = 1e-2
    mid_coils_edge = 1e-2
    outer_coils_edge = 1e-2
    trim_coils_edge = 2e-2
    
    center_coils_center = 0.5e-2 + center_coils_edge
    inner_coils_center = (center_coils_center + center_coils_edge + inner_coils_edge + 0.58e-2)
    mid_coils_center = (inner_coils_center + inner_coils_edge + mid_coils_edge + 1.683e-2)
    outer_coils_center = inner_coils_center + mid_coils_edge + outer_coils_edge + 18e-2
    trim_coils_center = outer_coils_center + outer_coils_edge + trim_coils_edge + 0.125e-2
    
    base_windings = [104848.32489944968]
    center_coil_windings =  base_windings[0]
    inner_coil_windings =   0.14015 * base_windings[0] / (2*inner_coils_edge)
    mid_coil_windings =  0.168 * base_windings[0] / (2 * mid_coils_edge)
    outer_coil_windings = 1.341e-3*base_windings[0] / (2*outer_coils_edge)
    trim_coil_windings = 3e-4 * base_windings[0] / (2*trim_coils_edge)

    #main magnet coils
    center_left_coil = Coil_form(inner_radius=center_coils_radii[0],
        outer_radius=center_coils_radii[1],
        left_edge=-center_coils_edge,
        right_edge=center_coils_edge,
        z_center=-center_coils_center,
        num_windings=center_coil_windings,
        current_per_wire=current_per_wire,
        name="left center coil")
        
    center_right_coil = Coil_form(inner_radius=center_coils_radii[0],
        outer_radius=center_coils_radii[1],
        left_edge=-center_coils_edge,
        right_edge=center_coils_edge,
        z_center=center_coils_center,
        num_windings=center_coil_windings,
        current_per_wire=current_per_wire,
        name="right center coil")
        
    inner_left_coil = Coil_form(inner_radius=inner_coils_radii[0],
        outer_radius=inner_coils_radii[1],
        left_edge=-inner_coils_edge,
        right_edge=inner_coils_edge,
        z_center=-inner_coils_center,
        num_windings=inner_coil_windings,
        current_per_wire=current_per_wire,
        name="left inner coil")

    inner_right_coil = Coil_form(inner_radius=inner_coils_radii[0],
        outer_radius=inner_coils_radii[1],
        left_edge=-inner_coils_edge,
        right_edge=inner_coils_edge,
        z_center=inner_coils_center,
        num_windings=inner_coil_windings,
        current_per_wire=current_per_wire,
        name="right inner coil")
        
    mid_left_coil = Coil_form(inner_radius=mid_coils_radii[0],
        outer_radius=mid_coils_radii[1],
        left_edge=-mid_coils_edge,
        right_edge=mid_coils_edge,
        z_center=-mid_coils_center,
        num_windings=mid_coil_windings,
        current_per_wire=-current_per_wire,
        name="left mid coil")

    mid_right_coil = Coil_form(inner_radius=mid_coils_radii[0],
        outer_radius=mid_coils_radii[1],
        left_edge=-mid_coils_edge,
        right_edge=mid_coils_edge,
        z_center=mid_coils_center,
        num_windings=mid_coil_windings,
        current_per_wire=-current_per_wire,
        name="right mid coil")
        
    outer_left_coil = Coil_form(inner_radius=outer_coils_radii[0],
        outer_radius=outer_coils_radii[1],
        left_edge=-outer_coils_edge,
        right_edge=outer_coils_edge,
        z_center=-outer_coils_center,
        num_windings=outer_coil_windings,
        current_per_wire=current_per_wire,
        name="left outer coil")

    outer_right_coil = Coil_form(inner_radius=outer_coils_radii[0],
        outer_radius=outer_coils_radii[1],
        left_edge=-outer_coils_edge,
        right_edge=outer_coils_edge,
        z_center=outer_coils_center,
        num_windings=outer_coil_windings,
        current_per_wire=current_per_wire,
        name="right outer coil")
        
    trim_coil_left = Coil_form(inner_radius=trim_coils_radii[0],
        outer_radius=trim_coils_radii[1],
        left_edge=-trim_coils_edge,
        right_edge=trim_coils_edge,
        z_center=-trim_coils_center,
        num_windings=trim_coil_windings,
        current_per_wire=-current_per_wire,
        name="left trim coil")
        
    trim_coil_right = Coil_form(inner_radius=trim_coils_radii[0],
        outer_radius=trim_coils_radii[1],
        left_edge=-trim_coils_edge,
        right_edge=trim_coils_edge,
        z_center=trim_coils_center,
        num_windings=trim_coil_windings,
        current_per_wire=-current_per_wire,
        name="right trim coil")

    list_coils= [center_left_coil,
        center_right_coil,
        inner_left_coil,
        inner_right_coil,
        mid_left_coil,
        mid_right_coil,
        outer_left_coil,
        outer_right_coil,
        trim_coil_left,
        trim_coil_right]

    field_profile = Field_profile(list_coils)

    return field_profile
    
def load_he6_coils(main_field,trap_strength=1e-3):

    #He6-CRES three-coil trap
    center_windings = 88 / (4.84e-3 * 2)
    edge_windings = 44 / (4.84e-3 * 2)
    current_per_wire = (0.24687194322084335 / 1e-3) * trap_strength * main_field
    

    center_coil = Coil_form(1.15e-2,1.35e-2,-4.84e-3,4.84e-3,0,center_windings,-current_per_wire,"center coil")

    left_edge_coil = Coil_form(1.15e-2,1.25e-2,-4.84e-3,4.84e-3,-4.5e-2,edge_windings,current_per_wire,"left coil")
    right_edge_coil = Coil_form(1.15e-2,1.25e-2,-4.84e-3,4.84e-3,4.5e-2,edge_windings,current_per_wire,"right coil")

    triple_field_profile = Field_profile([center_coil,
                                left_edge_coil,
                                right_edge_coil], trap_strength=1e-3)
                                
    return triple_field_profile
    
def load_he6_trap(main_field=1,trap_strength=1e-3):

    coil_list = load_he6_coils(main_field,trap_strength).get_coil_list()
    trap_profile = Trap_profile(coil_list,main_field)
    
    return trap_profile
