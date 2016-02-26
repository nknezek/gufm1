# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 11:52:39 2016

@author: nknezek
"""

from scipy.misc import factorial as _factorial
from scipy.special import lpmv as _lpmv
from numpy import sin as _sin
from numpy import cos as _cos

#%%
def Pml(x, l, m):
    """
    Associated Legendre Polynomial - Schmidt Quasi-Normalization
    ============================================================
    Returns the evaulated Associated Legendre Polynomial of degree n and order m at location x.
    
    This function evaluates the Associated Legendre Polynomials with Schmidt Quasi Normalization as defined in Schmidt (1917, p281).
    It uses the scipy built in associated legendre polynomials which have Ferrer's normalization and converts the normalization.    
    
    Inputs
    -------
    x:
        Location of evaluation
    l:
        Degree of associated legendre polynomial
    m:
        Order of associated legendre polynomial

    Returns 
    -------
    The value of the polynomial at location specified. (float)
    
    Associated Legendre Polynomial Normalizations:
    ------
    
    Schmidt Quasi-Normalized:
        P^m_l(x) = sqrt{2*(l-m)!/(l+m)!}(1-x^2)^{m/2}(d/dx)^2 P_l(x)
    
    Ferrer's (only for reference):
        P^m_n(x) = (-1)^m(1-x^2)^{m/2}(d/dx)^2 P_n(x)
        
    """
    return (2*_factorial(l-m)/_factorial(l+m))**0.5/(-1)**m*_lpmv(m,l,x) 


def Br_for_ml(r,th,ph,g,h,m,l, a=6371.2):
    """
    Calculates the Br contribution for one set of m,l, using the potential field.
        
    Inputs 
    ------
    r:
        radius location (km)
    th:
        latitude location (radians)
    ph:
        longitude location (radians)
    g:
        Gauss coefficient (cos term)
    h: 
        Gauss coefficient (sin term)
    m:
        Order of calculation
    l:
        Degree of calculation
    a: 
        Radius (km) at which Gauss coefficients are calculated

    Returns
    -------
    Br contribution in Tesla at a particular point from a particular degree and order.     
    """
    return (l+1.)*a**(l+2.)/abs(r)**(l+2.)*(g*_cos(m*ph) + h*_sin(m*ph))*Pml(_cos(th), l, m)
    
def Br(r,th,ph, g_dict, h_dict, l_max=None):
    '''
    Calculates the total radial magnetic field at a particular location, give a dictionary of gauss coefficients.
        
    Inputs 
    ------
    r:
        radius location (km)
    th: 
        latitude location (radians)
    ph:
        longitude location (radians)
    g_dict:
        dictionary of g (cos) Gauss coefficients, ordered as g[l][m].
    h_dict:
        dictionary of h (sin) Gauss coefficients, ordered as h[l][m]. h coefficients for m=0 should be explicitly included as 0.0
    l_max:
        maximum degree to use in calculation. By default uses all supplied degrees.
        
    Returns
    -------
    Total Br at a particular point (Tesla)
    '''
    if not l_max:
        l_max = len(g_dict)
    Br_sum = 0
    for l in range(1,l_max+1):
        for m in range(l+1):
            Br_sum += Br_for_ml(r,th,ph, g_dict[l][m], h_dict[l][m], m, l)
    return Br_sum

def read_gufm1_data(filename):
    '''
    Reads in gufm1 data from one timeknot and stores the data in two dictionaries for each (g,h) Gauss coefficient.
    
    Inputs
    ------
    filename:
        plain text file of gufm1 data, standard format as downloaded from website. Two header lines, interpreted as plain text.
    
    Returns
    -------
    g, h: 
        dictionaries of Gauss coefficients ordered as g[l][m] and h[l][m]
    
    '''
    with open(filename,'rb') as f:
        f.readline()
        l_max = int(f.readline().split()[0])
        data = []
        for line in f:
            for x in line.strip().split():            
                data.append(float(x))
        g = {}
        h = {}
        g[1] = {0:data[0]}
        g[1][1] = data[1]
        h[1] = {0:0, 1:data[2]}
        i = 3
        for l in range(2,l_max+1):
            g[l] = {}
            h[l] = {}
            g[l][0] = data[i]
            i += 1
            h[l][0] = 0.
            for m in range(1,l+1):
                g[l][m] = data[i]
                i += 1
                h[l][m] = data[i]
                i += 1
    return g, h