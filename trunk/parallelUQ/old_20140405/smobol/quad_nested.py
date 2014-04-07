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
import numpy as np


class QuadratureNested:
    """
    Nested quadrature rules in 1d, for use with sparse grids.  This is an
    abstract class containing common functionality; to define a particular rule
    define init_rule() (setting self.x and self.w, lists of nodes and weights
    respectively) and init_idx() (setting self.idx, the location of each node
    within the top-level rule).  See QuadraturePatterson and QuadratureCC for
    examples.
    """
    def __init__(self):
        self.init_rule()                # I.e. nodes and weights
        self.init_idx()
        self.init_dw()
        self.interpolate_pre()

    def integrate(self, level, f, fargs=()):
        """
        Most basic possible integration routine
          f     - a scalar function of one variable
          level - the level of rule, 0 to nlevel()-1
        """
        assert 0 <= level and level < self.nlevel()
        self.fdata = np.zeros(self.nnode(level))
        for i,x in enumerate(self.x[level]): self.fdata[i] = f(x, *fargs)
        return np.einsum('i,i', self.fdata, self.w[level])

    def interpolate_pre(self):
        """Pre-compute normalization factors for interpolation - store in cnorm"""
        self.cnorm = []
        for l in range(self.nlevel()):
            N = self.nnode(l)
            diff_mat  = np.array([self.x[l]])-np.array([self.x[l]]).T
            diff_mat += np.identity(N)  # Make the diagonal = 1
            self.cnorm.append(np.prod(diff_mat, axis=0))

    def interpolate(self, x, l, f, fargs=()):
        """
        x - scalar location [-1,1]
        l - level
        """
        self.fdata = np.array([f(xl, *fargs) for xl in self.x[l]])
        return np.einsum('i,i', self.fdata, self.lagrange(x,l))
            
    def lagrange(self, x, l):
        """Return Lagrange polys array (nnode(l),), level l, scalar location x.
        TODO: implement for vector x, to make interpolation more efficient."""
        N = self.nnode(l)
        c = np.array([x - self.x[l]]) * np.array([np.ones(N)]).T
        c = c * (1. - np.identity(N)) + np.identity(N)
        return np.prod(c, 1) / self.cnorm[l]
    
    def nlevel(self): return len(self.x)

    def nnode(self, level): return self.x[level].shape[0]
        
    def maxnnode(self): return self.x[-1].shape[0]

    def global_index_to_x(self, idx):
        """Convert scalar/list of global indices to node location(s)."""
        return self.x[-1][idx]

    def init_dw(self):
        """
        Init weights of rules \Delta_i(f) = U_i(f) - U_{i-1}(f).  Since nodes of
        U_{i-1} are a subset of those of U_i, the node locations of \Delta_i and
        U_i are identical.  Note U_{-1} = 0.
        """
        self.dw = [self.w[0]]
        for i in range(1, self.nlevel()):
            tmp = np.zeros(self.maxnnode())
            tmp[self.idx[i]] = self.w[i]
            tmp[self.idx[i-1]] -= self.w[i-1]
            self.dw.append( tmp[self.idx[i]] )
    
    def lagrange_delta(self, x, l):
        """
        Return differences between Lagrange polys of level l and level l-1,
        array (nnode(l),).  Location x.  These are basis functions for the
        component rules in the Smolyak algorithm.  Corresponds to init_dw(),
        which does the same for the weights.
        """
        if l == 0: return self.lagrange(x, 0)
        L = np.zeros(self.maxnnode())
        L[self.idx[l]]    = self.lagrange(x,l)
        L[self.idx[l-1]] -= self.lagrange(x,l-1)
        return L[self.idx[l]]

