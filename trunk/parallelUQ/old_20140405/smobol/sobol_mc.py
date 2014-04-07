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
import numpy as np


def sobol_variance_mc(f, dim, y, N=10000, monitor=False):
    """
    Compute Sobol variances D_y and D_y^{t_0t} of the function f of dim variables
    defined on the unit-hypercube.  List y is the index set - ordering is not
    important, for multiple index sets (e.g. all Sobol variances of 1 variable)
    multiple calls are required.  Method is Monte-Carlo with N samples.
    Return - (tuple) N        - Number of MC samples used
                     f0       - mean
                     D        - total variance 
                     D_y      - main effect variance of indices y
                     D_y_tot  - total effect variance of y
    Sensitivity indices are S_y        = D_y / D
                            S_y^{t_0t} = D_y^{tot} / D.
    Method described in:
    Sobol,"Global sensitivity indices for nonlinear mathematical models", 2001.
    """
    m = len(y)
                                        # Compliment of y
    z = sorted(list(set(range(dim)) - set(y)))
                                        # Return: N, f0, D, D_y, D_y_tot
    def postprocess(M, phi, phi2, psi, chi):
        f0 = phi / M
        return M, f0, phi2 / M - f0**2, psi / M - f0**2, chi / M
    if monitor:
        sys.stdout.write(('-'*100+'\n%10s %20s %20s %20s %20s\n'+'-'*100+'\n') %
                         ('N','f_0','D','D_y','D_y_tot'))
    phi, phi2, psi, chi = 0,0,0,0
    for i in xrange(1,N+1):
        eta   = np.random.random(m)     # Sample locations on [0,1]**dim
        etap  = np.random.random(m)
        zeta  = np.random.random(dim-m)
        zetap = np.random.random(dim-m)
        x1, x2, x3 = np.zeros(dim), np.zeros(dim), np.zeros(dim)
        x1[y] = eta;  x1[z] = zeta
        x2[y] = eta;  x2[z] = zetap
        x3[y] = etap; x3[z] = zeta
        phik = f(x1)                    # Eval f
        psik = phik * f(x2)
        chik = 0.5 * (phik - f(x3))**2
        phi  += phik                    # Monte-Carlo sum
        phi2 += phik**2
        psi  += psik
        chi  += chik
        if monitor and ((i > 9 and np.log10(i) % 1 < 1.e-8) or i == N):
            sys.stdout.write('%10d %20.8e %20.8e %20.8e %20.8e\n'
                             % postprocess(i, phi, phi2, psi, chi))
    return postprocess(N, phi, phi2, psi, chi)

#---------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------
def sobol_g_function(x, a):
    """
    Sobol's g-function (scalar)
      x - indep. variable, arbetrary dimension
      a - coefficients, same dimension as x and > 0
    """
    return np.prod((np.fabs(4.*x - 2.) + a) / (1. + a))

def sobol_g_function_exact(a):
    """Exact values of 1st-order Sobol indices (normalized) of g-function."""
    t = 1./(3*(1.+a)**2)
    return t / np.sum(t)


#---------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------
if __name__ == '__main__':
    dim = 6
    a = np.arange(1.,dim+1)
    def f(x): return sobol_g_function(x, a)
    
    D, D_y, D_y_tot = np.zeros(dim), np.zeros(dim), np.zeros(dim)
    for i in range(dim):
        print 'Dimension %d' % (i+1,)
        N, f0, D[i], D_y[i], D_y_tot[i] =sobol_index_mc(f, dim, [i], 
                                                        N=100000, monitor=True)

    exact = sobol_g_function_exact(a)
    print '\n%10s %14s %14s %14s' % ('Dim', 'Approx', 'Exact', 'Rel. error')
    for i in range(dim):
        S_y = D_y[i]/D[i]
        error = np.fabs(S_y - exact[i]) / exact[i]
        print '%10d %14.2e %14.2e %14.2e' % (i+1, S_y, exact[i], error)

