#   smobol - Sparse grid and Sobol approximations
#   This file is part of smobol.
#   Copyright 2013 Richard Dwight <richard.dwight@gmail.com>
#
#   smobol is free software: you can redistribute it and/or modify it under the
#   terms of the GNU Lesser General Public License as published by the Free
#   Software Foundation, either version 3 of the License, or (at your option)
#   any later version.
#
#   smobol is distributed in the hope that it will be useful, but WITHOUT ANY
#   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#   FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
#   more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with smobol.  If not, see <http://www.gnu.org/licenses/>.
import sys
from scipy.special import erfc
import numpy as np
try:
    import matplotlib.pyplot as plt
    plt_loaded = True
except ImportError:
    plt_loaded = False

import mylib


# Coefficients in Genz's test integrands.
c = []
w = []

#----------------------------------------------------------------------------
# N-dimensional test integrands.
#----------------------------------------------------------------------------
def f1(x):
    """Genz's "oscillatory" test integrand."""
    return np.cos(2.0 * np.pi * w[0] + np.sum(c * x))

def f2(x):
    """Genz's "product peak" test integrand."""
    return 1. / np.prod(1./c**2 + (x - w)**2)

def f3(x):
    """Genz's "corner peak" test integrand."""
    dim = len(x)
    return 1. / ((1. + np.sum(c * x))**(dim+1))

def f4(x):
    """Genz's "Gaussian" test integrand."""
    return np.exp( -np.sum((c * (x - w))**2) )

def f5(x):
    """Genz's "continuous" test integrand."""
    return np.exp( -np.sum(c * np.abs(x - w)) )

def f6(x):
    """Genz's "discontinuous" test integrand."""
    dim = len(x)
    if dim == 1:
        val = np.exp(where(x[0] <= w[0], np.sum(c * x), 0.0))
    else:
        if(x[0] <= w[0] and x[1] <= w[1]): val = np.exp(np.sum(c * x))
        else: val = 0.
    return val

#----------------------------------------------------------------------------
# Corresponding exact solutions.
#----------------------------------------------------------------------------
def f1_exact(dim):
    arg = np.sum(c)
    prod = reduce(lambda x, y: x*y, np.sin(0.5 * c) / c)
    return 2.0**dim * np.cos(2.0 * np.pi * w[0] + 0.5 * arg) * prod;

def f2_exact(dim):
    return np.prod(c * (np.arctan(c * (1. - w)) + np.arctan(c * w)))

def f3_exact(dim):
    raise NotImplementedError() # TODO

def f4_exact(dim):
    return np.prod(np.sqrt(np.pi)/(2*c) * (np.erfc(-c * w) - np.erfc(c * (1.-w))))

def f5_exact(dim):
    return np.prod(1./c * (2. - np.exp(-c*w) - np.exp(c*(w-1.))))

def f6_exact(dim):
    if dim == 1:
        val = np.prod((np.exp(c * w) - 1.) / c)
    else:
        c1 = (np.exp(c * w) -1.) / c
        c2 = (np.exp(c) -1.) / c
        val = np.prod(c1[0:2]) * np.prod(c2[2:])
    return val


#----------------------------------------------------------------------------
# 
#----------------------------------------------------------------------------
def get_fn(fnnum):
    """Select a test function by index."""
    return eval('f'+str(fnnum))

def get_fn_exact(fnnum):
    """Select a test function exact solution by index."""
    return eval('f'+str(fnnum)+'_exact')

def plot_testfn():
    """Make contour plots of all 6 test functions in 2d."""
    assert plt_loaded, 'Matplotlib not installed'
    dim = 2
    global c, w
    c = np.array([0.5] * dim); c = c / sum(c) * 9.
    w = np.array([0.5] * dim)

    xi = np.linspace(0., 1., 100)
    xx = mylib.meshgrid_flatten(xi, xi)
    
    fig = plt.figure(figsize=(12, 8))
    for i in range(1, 7):
        fn = get_fn(i)
        ax = fig.add_subplot(2, 3, i)
        F = np.zeros(xx.shape[0])
        for i, x in enumerate(xx):
            F[i] = fn(np.array(x))
        F.reshape(xi.size, xi.size)
        ax.contour(xi, xi, F.reshape(xi.size, xi.size).T)

    fig.savefig('test_genz.contours.pdf')


#----------------------------------------------------------------------------
# 
#----------------------------------------------------------------------------
if __name__ == '__main__':
    plot_testfn()

