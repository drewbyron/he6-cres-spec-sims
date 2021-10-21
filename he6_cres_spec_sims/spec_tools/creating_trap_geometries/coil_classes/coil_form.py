import math
import numpy as np
from scipy.integrate import romberg
from scipy.misc import derivative


class Coil_form:
    """
    Generates object that contains geometric data of current coil. Note that num_windings is the number of windings per unit length (winding density).
    """
    
    def __init__(self,inner_radius,outer_radius,left_edge,right_edge,z_center,num_windings,current_per_wire,name=None):
    
        if not name == None:
            self._name = str(name)
        else:
            self._name = name
            
        self._inner_radius = inner_radius
        self._outer_radius = outer_radius
        self._left_edge = left_edge
        self._right_edge = right_edge
        self._z_center = z_center
        self._num_windings = num_windings
        self._current_per_wire = current_per_wire
        
        while self._inner_radius >= self._outer_radius:
        
            print("ERROR: inner radius greater than outer radius")
            self._inner_radius = input("Reset inner radius: ")
        
    def get_coil_dimensions(self):
        """
        Returns coil dimensions.
        """
        
        return self._inner_radius,self._outer_radius,self._left_edge,self._right_edge
        
    def get_coil_properties(self):
        """
        Returns coil center z-coordinate, winding length density, and current per wire.
        """
        
        return self._z_center,self._num_windings, self._current_per_wire
                
    def set_coil_dimensions(self,inner_radius,outer_radius,left_edge,right_edge):
        """
        Set coil dimensions.
        """
        self._inner_radius = inner_radius
        self._outer_radius = outer_radius
        self._left_edge = left_edge
        self._right_edge = right_edge
        
    def set_coil_properties(self,z_center,num_windings,current_per_wire):
        """
        Set coil center z-coordinate, winding length density, and current per wire.
        """
        self._z_center = z_center
        self._num_windings = num_windings
        self._current_per_wire = current_per_wire
        
    def field_values(self,position,sheets=20,field_coordinates="Cartesian",calculate_Br=True):
        """
        Calculates the magnetic field T in coordinates field_coordinates at position = (radius,phi,z) in cylindrical coordinates
        with phi in radians. Can choose field_coordinates 'Cartesian' or 'cylindrical'.
        """
        
        mu0 = 4 * math.pi * 1e-7
        current_density = self._num_windings * self._current_per_wire
        
        #break down position vector
        radius = position[0]
        phi = position[1]
        zpos = position[2]
        
        #shift coordinate system for integral calculation
        zpos = zpos - self._z_center
        
        if radius == 0:
            
            def logterm(z,zedge,Ro,Ri):
                    numerator = Ro + math.sqrt(Ro**2 + (z-zedge)**2)
                    denominator = Ri + math.sqrt(Ri**2 + (z-zedge)**2)
                    value = (z-zedge)*math.log(numerator / denominator)
       
                    return value
       
            Br = 0
            Bz = (mu0 / 2)* (current_density /(self._outer_radius - self._inner_radius))*(logterm(zpos,self._left_edge,self._outer_radius,self._inner_radius)-logterm(zpos,self._right_edge,self._outer_radius,self._inner_radius))
            
        else:
            
            def Q_func(z, zedge, r, R, theta):
                return math.sqrt((z-zedge)**2 + r**2 + R**2 - 2*r*R*math.cos(theta))
                
                
            def Br_integrand(z,z1,z2,r,R,theta):
            
                try:
                    Q_one = Q_func(z,z1,r,R,theta)
                    Q_two =  Q_func(z,z2,r,R,theta)
                    
                    first_term_numerator = R*math.cos(theta) * (1 + (z-z2)/Q_two)
                    first_term_denominator = (z-z2) + Q_two
                    first_term = first_term_numerator/first_term_denominator
                    
                    second_term_numerator = -R*math.cos(theta) * (1 + (z-z1)/Q_one)
                    second_term_denominator = (z-z1) + Q_one
                    second_term = second_term_numerator/second_term_denominator
                    
                    integrand = first_term + second_term
                    
                except:
                    print("WARNING: Unhandled value for Br integrand at radius = {}, theta = {}, z = {}".format(r,theta,z))
                    print("Returning 0 for integrand; may affect accuracy")
                    
                    
                    integrand = 0
                
                return integrand
                
            def Bz_integrand(z,z1,z2,r,R,theta):
            
                try:
    
                    Q_one = Q_func(z,z1,r,R,theta)
                    Q_two =  Q_func(z,z2,r,R,theta)
                    
                    first_term_numerator = R*math.cos(theta) * (2 * r - 2 * R * math.cos(theta))
                    first_term_denominator = 2 * Q_one * (z-z1 + Q_one)
                    first_term = first_term_numerator / first_term_denominator
    
    
                    second_term_numerator = -R*math.cos(theta) * (2 * r - 2 * R * math.cos(theta))
                    second_term_denominator = 2 * Q_two * (z-z2 + Q_two)
                    second_term = second_term_numerator / second_term_denominator
    
                    third_term = R * math.cos(theta) / r * math.log((z-z1+Q_one)/(z-z2+Q_two))
    
                    integrand = first_term +  second_term + third_term
                    
                except:
                    print("WARNING: Unhandled value for Bz integrand at radius = {}, theta = {}, z = {}".format(r,theta,z))
                    print("Returning 0 for integrand; may affect accuracy")
                    integrand = 0
    
                return integrand
            
            current_per_sheet = current_density / sheets
            Br = 0
            Bz = 0
            curr_sheets = 0
            for k in range(sheets):
            
                def Br_field_function(theta):
                    return Br_integrand(zpos,self._left_edge,self._right_edge,radius,self._inner_radius + k * (self._outer_radius - self._inner_radius)/(sheets-1),theta)
    
                def Bz_field_function(theta):
                    return Bz_integrand(zpos,self._left_edge,self._right_edge,radius,self._inner_radius + k * (self._outer_radius - self._inner_radius)/(sheets-1),theta)
                    
                if calculate_Br == True:
                    Br = Br + (mu0 / (2*math.pi)) * current_per_sheet * romberg(Br_field_function,0,math.pi,divmax=20)
                else:
                    Br = 0
                
                Bz = Bz + (mu0 / (2*math.pi)) * current_per_sheet * romberg(Bz_field_function,0,math.pi,divmax=20)
              
        if field_coordinates == "cylindrical":
            return (Br, 0 , Bz)
        
        elif field_coordinates == "Cartesian":
              
            Bx = Br * math.cos(phi)
            By = Br * math.sin(phi)
    
            return (Bx,By,Bz)
            
        else:
            print("ERROR: {} not a valid coordinate system".format(return_coordinates))
            
    def field_strength(self,radius,zpos):
        """
        Calculates Bz magnetic field strength in T at radius, zpos.
        """
    
        position = (radius,0,zpos)
        total_field = self.field_values(position,calculate_Br=False)[2]
     
        return total_field
        
    def field_derivative(self,radius,zpos,deriv_order=1,dx=1e-6):
        """
        Calculates dBz / dz derivative on z-axis.
        """

        def f(x):
            return self.field_strength(radius,x)

        return derivative(f,zpos,dx=dx,n=deriv_order,order=(deriv_order*2 + 1))
   
    @property
    def inner_radius(self):
        return self._inner_radius
    
    @property
    def outer_radius(self):
        return self._outer_radius
           
    @property
    def left_edge(self):
        return self._left_edge
           
    @property
    def right_edge(self):
        return self._right_edge
           
    @property
    def z_center(self):
        return self._z_center
        
    @z_center.setter
    def z_center(self,value):
        self._z_center = value
        
    @property
    def num_windings(self):
        return self._num_windings
        
    @num_windings.setter
    def num_windings(self,value):
        self._num_windings = value
        
    @property
    def current_per_wire(self):
        return self._current_per_wire
        
    @current_per_wire.setter
    def current_per_wire(self,value):
        self._current_per_wire = value
        
    @property
    def name(self):
        return self._name
        
    @name.setter
    def name(self,value):
        self._name = str(value)

