from spec_tools.spec_gen import spec_gen as sg
from spec_tools.spec_calc import spec_calc as sc

class He6_energy_spec:
    """
    Generates a normalized decay spectrum for He6 with specified step
    size in eV and choice of Fierz coefficient b.
    """
    def __init__(self,step_size=100,b=0):
      
        self.xvals=[]
        self.yvals=[]
      
        spec_range = 3.508e6 - 0.0024
        num_bins = int(spec_range / step_size)
        
        normalization = sg.spectrum_prob_norm(b)
        
        for k in range(num_bins-1):
      
            self.xvals.append(0.0024 + k * step_size)
            self.yvals.append(sg.spectrum_prob_range(self.xvals[k],self.xvals[k] + step_size,b,normalization))
            
            k = k + 1
      
        self._size = num_bins
        self._step = step_size
        self._b = b
        
    def energy_vals(self):
        """
        Returns (x_values,y_values) where x_values are in electron volts (eV).
        """
    
        xval_energy = [x for x in self.xvals]
        yval_energy = [y for y in self.yvals]
        
        return [xval_energy, yval_energy]
        
    def freq_vals(self,field,minFreq=0,maxFreq=0):
        """
        Returns (x_values, y_values) where x_values converted to frequency.
        """
        x_freq_vals = [sc.energy_to_freq(energy,field) for energy in self.xvals[1:]]
        y_freq_vals = [y for y in self.yvals[0:-1]]
        
        if not (minFreq == 0 and maxFreq == 0):
            xypairs = [pair for pair in zip(x_freq_vals,y_freq_vals) if pair[0] >= minFreq and pair[0] <= maxFreq]
            x_freq_vals = [pair[0] for pair in xypairs]
            y_freq_vals = [pair[1] for pair in xypairs]
            
        x_freq_vals.reverse()
        y_freq_vals.reverse()
        
        return [x_freq_vals,y_freq_vals]
    
    @property
    def size(self):
        return self._size
      
    @property
    def step(self):
        return self._step
        
    @property
    def b(self):
        return self._b
