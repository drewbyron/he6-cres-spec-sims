import math
import numpy as np
from scipy.misc import derivative
from scipy.optimize import fmin

from he6_cres_spec_sims.spec_tools.coil_classes.field_profile import Field_profile

class Trap_profile(Field_profile):
    """
    Generates object that calculates field strength of magnetic mirror trap of given coils in list_coils. Returns None if list_coils does not represent a proper trap.
    """
    
    #elements in Field_profile
    #self._list_coils
    #self._main_field
    #get_coil_list(self)
    #set_coil_list(self,list_coils)
    #field_values(self,position,field_coordinates="Cartesian")
    #field_grad(self,position,deriv_order=1,dx=1e-6,grad_coordinates="Cartesian")
    #field_strength(self,position,field_coordinates="Cartesian")
    #field_derivative(self,position,deriv_order=1,dx=1e-6)
    #@property main_field(self), main_field(self,value)
    
    
    def __init__(self,list_coils,main_field=0,field_scales=True):
    
        
        super().__init__(list_coils,main_field)
        self._field_scales = field_scales
        
        if not main_field > 0:
            print("WARNING: Main field not greater than 0")
            print("Not a valid trap...")
            self._is_trap = False
            return None
            
        elif not (self.field_derivative(0,0) == 0
                and self.field_derivative(0,0,2) > 0):
            print("WARNING: Given field profile does not have a local minimum at z=0")
            print("Not a valid trap...")
            self._is_trap = False
            return None
            
        self._is_trap = True
        self._trap_width = self.trap_width_calc()
        
    def trap_width_calc(self):
        """
        Calculates the trap width of the object trap_profile.
        """
    
        field_func = self.field_strength
        def func(z):
            return -1 * field_func(0,z)
        
        maximum = fmin(func,0,xtol=1e-12)[0]
        print("Trap width: ({},{})".format(-maximum,maximum))
        print("Maximum Field: {}".format(-1 * func(maximum)))
    
        trap_width = (-maximum,maximum)
        return trap_width
        

    @property
    def main_field(self):
        return self._main_field
        
    @main_field.setter
    def main_field(self,value):
        
        if not value > 0:
            print("ERROR: New main field must be greater than 0")
            return
        
        else:
            if self._field_scales == True:
                for coil in self._list_coils:
                    base_current = coil.current_per_wire / self._main_field
                    coil.current_per_wire = base_current * value
            self._main_field = value
        
    @property
    def is_trap(self):
        return self._is_trap
        
    @property
    def field_scales(self):
        return self._field_scales
        
    @property
    def trap_width(self):
        return self._trap_width
        
    
