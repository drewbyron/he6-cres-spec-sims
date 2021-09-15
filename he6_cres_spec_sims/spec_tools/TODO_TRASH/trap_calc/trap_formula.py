import numpy as np
from numpy.linalg import solve as ln_solve
from numpy.linalg import inv as ln_inv

from scipy.optimize import curve_fit
from spec_tools.spec_calc import spec_calc as sc

    
class Trap_formula():
    """
    Calculates trapping formula for he6 betas for a given trap radius
    as function of cyclotron frequency and main field strength.
    """
    def __init__(self, xdata, ydata, main_field, trap_radius):
    
        self.xdata = xdata
        self.ydata = ydata
        self._main_field = main_field
        self._trap_radius = trap_radius
        
        fit_xdata = [sc.cyc_radius(sc.freq_to_energy(xval,main_field),main_field,90)/trap_radius for xval in xdata]
        fit_ydata = [vals for vals in ydata]
        
        self._parameters = formula_fit(fit_xdata,fit_ydata)
        
    
    def trap_prob(self,freq,main_field,trap_radius):
        """
        Calculates the trapping probability of a He6 decay beta given
        freq in Hz, main_field in T, and trap_radius in meters.
        """
        a,b,c,d = tuple(self._parameters)
        cyc_radius = sc.cyc_radius(sc.freq_to_energy(freq,main_field),main_field,90)
        
        if cyc_radius == 0:
            print("WARNING: electron likely has cyclotron frequency greater than maximum possible; trapping probability set to zero")
            trap_probability = 0
        elif cyc_radius < trap_radius:
            radii_ratio = cyc_radius / trap_radius
            trap_probability = a/radii_ratio + b*radii_ratio + c*np.power(radii_ratio,2) + d
        else:
            print("WARNING: current cyclotron radius exceeds trap radius")
            print("Trapping probability set to zero")
            trap_probability = 0
     
        return trap_probability
        
    def print_formula(self):
        a,b,c,d = tuple(self._parameters)
        formula = '{:.4e}/x + {:.4e}x + {:.4e}x^2 + {:.4e}'.format(a,b,c,d)
        return formula
        
    @property
    def main_field(self):
        return self._main_field
        
    @property
    def trap_radius(self):
        return self._trap_radius
    
    @property
    def parameters(self):
        return self._parameters
        
#non-class functions used to find fit for trapping formula
def model_func(x,a,b,c,d):
    """
    model function for trap probability fitting
    """
    return a/x + b*x + c*np.power(x,2) + d

def formula_fit(xdata,ydata):
    """
    Calculates the paramters for model_func that fits the data
    """

    popt, pcov = curve_fit(model_func,xdata,ydata)
    print("covariance matrix: ",pcov)
    print("parameters:",popt)
    return popt
    
    
        
        
        
        
        
        
    

    
