import numpy as np
import math
import cmath

def gmma(y):

    z = 0
    gmma_ = 0

    f = 1.0
    n = 0
    c = [0.99999999999999044, 0.42278433510233479,
    0.41184033016678129, 0.08157692612415546,
    0.07424891541944474, -0.00026618659495306,
    0.01114971433577893, -0.00283646253037282,
    0.00206109185022554, -0.00083756468513517,
    0.00037536505226307, -0.00012141734870632,
    0.00002798328899383, -0.00000303019081028]
    
    while (y + n < 2 or y + n >= 3):
    
        if (y + n < 2):
            f = f * (y + n)
            n = n + 1
        if (y + n >= 3):
            n = n - 1
            f = f * (y + n)
    
    if (y + n == 2):
        gmma_ = 1.0
    else:
        gmma_ = c[13]
        z = y + n - 2
        for ii in range(13):
            i = 12 - ii
            gmma_ = gmma_ * z + c[i]
    
    if ( n < 0):
        gmma_ = gmma_ * f
    
    if ( n > 0):
        if (abs(f) < 1.0e-30):
            gmma_ = 1.0e30
        else:
            gmma_ = gmma_ / f
            
    return gmma_
    
def clngam(x):

    y = x
    g = 1.0 + 0.0j
    z = 0.0 + 0.0j
    clngam_ = 0.0 + 0.0j
    
    while (y.real <= 7):
        g = g * y
        y = y + 1.0
        
    z= -y**(-2)
    c=.398942280401433e0
    d=.295506535947712e-1
    a1=.282001658833287e1
    a2=.940005529444291e-1
    a3=.268573008412655e-1
    a4=.201429756309491e-1
    a5=.284850160437664e-1
    a6=.648894925920085e-1
    a7=.216924352948683e0
    clngam_ = (y - 0.5) * np.log(y) - y - np.log(c * g) + (a1 + z * (a2 + z * (a3 + z * (a4 + z * (a5 + z * (a6 + z * (a7 + z))))))) * d / y

    return clngam_
    
def coulomb_function(z,e,a):
    """
    Coulomb function for He6 beta decay spectrum calculation.
    """
    zi = 0.0 + 1.0j
    clg = 0. + 0.j
    
    xme = 0.510999
    alpha = 0.007297203
    
    e = xme * e
    xke = (e ** 2 - xme ** 2) ** (0.5)
    rad = (1.2 / 197.3) * (a ** (1/3))
    gam1 = (1.0 - (alpha * z) ** 2) ** (0.5)
    eta = z * alpha * e / xke
    
    clg = clngam(gam1 + zi * eta)
    z1_comp = cmath.exp(clg)
    z1= abs(z1_comp)
    z2= gmma(2 * gam1 + 1.)
    rat= z1 / z2
    expo= math.exp(3.1415926 * eta)

    return rat * rat * 4. * pow(( 2. * xke * rad), 2. * gam1 - 2.) * expo
