# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 11:52:39 2016

@author: nknezek
"""
import numpy as np
import scipy as sp
from scipy.misc import factorial
from scipy.special import lpmv
from numpy import sin, cos
import matplotlib.pyplot as plt
import cPickle as pkl

#%%
def Pml(x, l, m):
    """
    **Associated Legendre Polynomial - Schmidt Quasi-Normalization**
    
    This function evaluates the Associated Legendre Polynomials with Schmidt Quasi Normalization as defined in Schmidt (1917, p281).
    It uses the scipy built in associated legendre polynomials which have Ferrer's normalization and converts the normalization.    
    
    Returns the evaulated polynomial of degree n and order m at location x.

    x:
        Location of evaluation
    l:
        Degree of associated legendre polynomial
    m:
        Order of associated legendre polynomial
        
    **Associated Legendre Polynomial Normalizations:**
    
    Schmidt Quasi-Normalized:
        P^m_l(x) = sqrt{2*(l-m)!/(l+m)!}(1-x^2)^{m/2}(d/dx)^2 P_l(x)
    
    Ferrer's (only for reference):
        P^m_n(x) = (-1)^m(1-x^2)^{m/2}(d/dx)^2 P_n(x)
        
    """
    return (2*factorial(l-m)/factorial(l+m))**0.5/(-1)**m*lpmv(m,l,x) 


def Br_for_ml(r,th,ph,g,h,m,l, a=6371.2):
    """
    Potential Field
    """
    return (l+1.)*a**(l+2.)/abs(r)**(l+2.)*(g*cos(m*ph) + h*sin(m*ph))*Pml(cos(th), l, m)
    
def Br(r,th,ph, g_dict, h_dict, l_max=None):
    if not l_max:
        l_max = len(g)
    Br_sum = 0
    for l in range(1,l_max+1):
        for m in range(l+1):
            Br_sum += Br_for_ml(r,th,ph, g[l][m], h[l][m], m, l)
    return Br_sum

def read_gufm1_data(filename):
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
#%%

filename  = 'gufm1_1990.txt'
g,h = read_gufm1_data(filename)
#%%
Nth = 200
Nph = 180

th = np.linspace(0,np.pi,Nth)
ph = np.linspace(0, 2.*np.pi*(Nph-1)/Nph,Nph)
Brth = []
for t in th:
    Brph = []
    for p in ph:
        Brph.append(Br(3480, t, p, g, h, l_max=14)/1e6)
    Brth.append(np.average(np.abs(Brph)))
pkl.dump(Brth,open('Br_long_avg200.p','wb'))

#%%
plt.plot(th*180/np.pi, np.abs(Brth))
plt.grid()
plt.title('GUFM1 Br magnitude at Core surface, Longitude-averaged')
plt.xlabel('colatitude (degrees)')
plt.ylabel('Br magnitude (mT)')
plt.ylim([0,0.7])
#%%
g_avg = []
h_avg = []
for i in range(1,len(g)+1):
    gtmp = []
    htmp = []
    for j in range(i):
        gtmp.append(np.abs(g[i][j]))
        htmp.append(np.abs(h[i][j]))
    g_avg.append(np.average(gtmp))
    h_avg.append(np.average(htmp))
plt.semilogy(g_avg)
plt.semilogy(h_avg)
plt.title('GUFM1 Power Spectrum')
plt.xlabel('Degree (l)')
plt.ylabel('Average magnitude of coefficients')
plt.legend(['g','h'], loc='best')
plt.grid()