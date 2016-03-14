# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 11:52:39 2016

@author: nknezek
"""

from scipy.misc import factorial as _factorial
from scipy.special import lpmv as _lpmv
from numpy import sin as _sin
from numpy import cos as _cos
from numpy import zeros as _zeros
import os

data_file = os.path.dirname(os.path.abspath(__file__)) + '/data/gufm1.txt'

def read_gufm1tk_data(filename):
    '''
    Reads in gufm1 data from one timeknot and stores the data in two dictionaries for each (g,h) Gauss coefficient.

    Inputs
    ------
    filename:
        plain text file of gufm1 data, standard format as downloaded from website. Two header lines, interpreted as plain text.

    Returns
    -------
    l_max:
        spherical harmonic degree of model (14)
    data:
        gauss coefficients in raw ordering at timeknot
    '''
    with open(filename,'rb') as f:
        f.readline()
        l_max = int(f.readline().split()[0])
        data = []
        for line in f:
            for x in line.strip().split():
                data.append(float(x))
    return l_max, data

def read_gufm1tk_to_gh(filename):
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
    l_max, data = read_gufm1tk_data(filename)
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

def read_gufm_all(filename=data_file):
    '''

    Parameters
    ----------
    filename

    Returns
    -------

    '''
    with open(filename,'rb') as f:
        f.readline()
        line1 = f.readline().split()

        l_max = int(line1[0])
        nspl = int(line1[1])
        n = l_max*(l_max+2)

        gt = _zeros(n*nspl)
        tknts = _zeros(nspl+4)
        tknts[:3] = [float(x) for x in line1[2:]]
        ti = 3
        gi = 0
        for line in f:
            if ti+4 <= len(tknts):
                tknts[ti:ti+4] = [float(x) for x in line.split()]
                ti += 4
            else:
                gt[gi:gi+4] = [float(x) for x in line.split()]
                gi += 4
    gt_out = gt.reshape(n, nspl, order='F')
    return gt_out, tknts, l_max, nspl

def interval(tknts, time):
    '''
    Calculates nleft: the index of the timeknot on the left of the interval
        tknts[nleft] < tknts[nleft+1]
        tknts[nleft] <= time <= tknts[nleft+1]

    Parameters
    ----------
    tknts:
        a numpy array containing the timestamps for all knots in the model
    time:
        the time to calculate the field

    Returns
    -------
    the index of the time knot on the left of the interval
    '''
    if (time >= tknts[3] and time <= tknts[-4]):
        for n in range(3,len(tknts)):
            if time >= tknts[n]:
                nleft = n
            else:
                break
    else:
        raise IndexError("The time you've chosen is outside this model")
    return nleft

def bspline(time, tknts, jorder=4):
    '''
    Calculates B-spline and time knot index location for time t.

    Parameters
    ----------
    time:
        time to calculate
    tknts:
        array of time-knots
    jorder:
        order of b-splines

    Returns
    -------
    nleft:
        index of the time knot on the left of the interval (tknts[nleft] <= time <= tknts[nleft+1])
    spl:
        array of dimension jorder (default 4) containing the spline factors at time t.
    '''

    nleft = interval(tknts, time)

    deltal = _zeros(jorder-1)
    deltar = _zeros(jorder-1)
    spline = _zeros(jorder)

    spline[0] = 1.0
    for j in range(jorder-1):
        deltar[j] = tknts[nleft+j+1] - time
        deltal[j] = time - tknts[nleft-j]
        saved = 0.0
        for i in range(j+1):
            term = spline[i]/(deltar[i]+deltal[j-i])
            spline[i] = saved + deltar[i]*term
            saved = deltal[j-i]*term
        spline[j+1] = saved
    return nleft, spline

def calculate_gt_raw(gt, spl, nleft, l_max=14, jorder=4):
    '''
    Calculates the Gauss Coefficients in raw ordering given the parameters calculated by inverval() and bspline().

    Parameters
    ----------
    gt:
        raw data from gufm1 (n x nspl numpy array)
    spl:
        B-spline basis (jorder numpy array)
    nleft:
        coordinate of the timeknot to the left of desired time
    l_max:
        spherical harmonic degree included in model (14)
    jorder:
        order of B-splines (4)
    Returns
    -------
        Gauss Coefficients for time in raw ordering.
    '''
    n = l_max*(l_max+2)
    g_raw = _zeros(n)
    for k in range(n):
        for j in range(jorder):
            g_raw[k] += spl[j]*gt[k,j+nleft-4]
    return g_raw

def convert_data_to_gh(data, l_max=14):
    '''
    Converts data computed for a time to g, h dictionaries

    Inputs
    ------
    data:
        numpy array of data, standard ordering as on single-time data files from website.
    l_max:
        spherical harmonic degree included in model (14)
    Returns
    -------
    g, h:
        dictionaries of Gauss coefficients ordered as g[l][m] and h[l][m]

    '''
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

def get_gh_at_t(time, filename=data_file,  jorder=4):
    gt, tknts, l_max, nspl = read_gufm_all(filename)
    nleft, spl = bspline(time, tknts, jorder=jorder)
    data = calculate_gt_raw(gt, spl, nleft, l_max=l_max, jorder=jorder)
    g_dict ,h_dict = convert_data_to_gh(data)
    return g_dict, h_dict

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

def Rl(l, g, h, r=6371.2, a=6371.2):
    '''
    Calculates the mean-square field for a particular degree (l)
    '''
    Rsum = 0
    for m in g[l].iterkeys():
        Rsum += (l+1)*(g[l][m]**2+h[l][m]**2)
    return Rsum*(a/r)**(2.*l+4.)

def Rl_list(g,h,r=6371.2, a=6371.2):
    '''
    Calculates the mean-square field for all degrees (l)
    '''
    Rll = []
    for l in g.iterkeys():
        Rll.append(Rl(l, g, h, r=r, a=a))
    return Rll

def Br_rms_sq(Rl_list):
    '''
    Calculates the means-square radial field at a particular radius
    '''
    Br_sum = 0.
    for l,Rl in zip(range(1,len(Rl_list)+1), Rl_list):
        Br_sum += (l+1.)/(2.*l+1)*Rl
    return Br_sum

