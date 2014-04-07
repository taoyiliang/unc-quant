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
import sys, copy, itertools

import numpy as np
try:
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    plt_loaded = True
except ImportError:
    plt_loaded = False

import mylib
from quad_patterson import QuadraturePatterson
from quad_cc import QuadratureCC

class LevelSet:
    """
    Level-set describing a general sparse grid.  Each "level" (array (dim,),
    dtype=int8) in the set corresponds to one difference rule in the Smolyak
    sparse-grid sum (i.e. one tensor-product of 1d rules (\Delta_i(f) = U_i(f) -
    U_{i-1}(f)).  Standard Smolyak sparse-grid can be described by setting-up
    the level-set s.t. it contains all levels k with |k| < l+dim-1.  This is
    achieved with init_simplex(l).  Level-sets are "admissible" (i.e. result in
    valid quadrature rules) iff they are 

    Adaptivity implementation following: Gerstner + Griebel, Dimension-Adaptive
    Tensor-Product Quadrature, Computing 71(1), 2003 (note, there called "index
    set").  New levels are added to the level-set based on an error-estimation,
    with attention to preserving the admissibility of the set.  Pays lots of
    attention to memory and complexity for large dim.  Notation follows that
    paper too.

    Definitions:
      - Level (vector)     = array (dim,) describing quadrature level in each
                             dimension, minimum (non-zero) level == 1.  A
                             sparse-grid rule is constructed from a sum of these
                             rules.
      - Level set I        = all levels contained in the sparse-grid rule

      The following are only needed for adaptivity:
      - I-idx              = a (scalar) index, a location in I
      - Active index set A = (a set of I-idx) levels in I whose integrals have
                             been computed, but all of whose forward neighbours
                             have not been computed, "vanguard" indices
      - Old index set O    = (a set of I-idx) all other levels in I, i.e I \ A
                             The active/old distinction is needed only for 
                             adaptivity, not relevant for description of rule
    """
    def __init__(self, dim, level, fill='simplex'):
        """Initialize level-set according to fill-pattern and level."""
        self.dim=dim
                                        # Main info, list of all indexes in
                                        # the set.  We then refer to them by
                                        # reference.  This has to grow, so
                                        # it's a list rather than ndarray
        self.I  = [np.ones(self.dim,dtype=np.int8)]
        self.A  = set([0])              # Active multi-indices (on the
                                        # forward bdry) described by their
                                        # location in I
        self.O  = set()                 # Old multi-idxs (by location in I)
        if   fill=='simplex':     self.init_simplex(level)
        elif fill=='full_factor': self.init_full_factor(level)
        elif fill=='one_factor':  self.init_one_factor(level)
        elif fill=='none':        pass
        else: raise ValueError('Unknown fill pattern: %s' % fill)
        # TODO: adaptivity
        #self.Nf = []                   # Forward and backward neighbours
        #self.Nb = []                   # (location in I)
        #self.G  = []                   # Error estimates

    def init_one_factor(self, l):
        """Initialize with one-factor analysis on each variable."""
        self.active_to_old(0)
        for i in range(self.dim):
            for j in range(2, l+1):
                k = np.ones(self.dim); k[i] = j
                self.add(k)

    def init_simplex(self, l):
        """
        Initialize with the "standard" sparse-grid index set of level l>0: 
        |k| <= l+dim-1.  This choice has some nice theoretical error properties.
        Recursive implementation.
        """
        assert l > 0, 'Level must be >0: l = %d' % l
        def r(m, d):
            if d == 1: return np.array([np.arange(m+1)], dtype=np.int8).T
            for n in range(m+1):
                A = r(m-n, d-1)
                b = np.ones((A.shape[0], 1), dtype=np.int8)*n
                bA = np.hstack((b, A))
                if n == 0: out = bA
                else:      out = np.vstack((out, bA))
            return out
        self.I = r(l-1, self.dim) + 1   # Directly set multi-indices
        self.A = set(); self.O = set()  # Directly set active/old
        for i, midx in enumerate(self.I):
            if np.sum(midx) == l+self.dim-1: self.A.add(i)
            else:                            self.O.add(i)

    def init_full_factor(self, l):
        """
        Initialize with the full tensor-product (non-sparse) of level l>0.  This
        is just for reference purposes, this call will be way too expensive in
        high-dimensions, and the implementation of this rule via sparse-grids is
        ridiculously convoluted.
        """
        self.active_to_old(0)
        Itmp = mylib.meshgrid_flatten( *(range(1,l+1),)*self.dim )
        for k in Itmp[1:]:              # 1st entry is already there
            i = self.add(k)
            if max(k) < l: self.active_to_old(i)

    def reduce_axes(self, axes):
        """
        Create new LevelSet by reducing this one along specified dimensions.
        Corresponds to equivalent, lower- dimensional sparse-grid rule.  This
        LevelSet remains untouched.
          axes - dimensions to remove from the levelset
        Return:
          ls_r - LevelSet of reduced dimension
        """
        axes_r = mylib.complement(axes, self.dim)
        dim_r = len(axes_r)
        ls_r = LevelSet(dim_r, 0, fill='none')
                                        # Delete unwanted dimensions and remove
                                        # duplicate multi-indices (lose order)
        ls_r.I = np.array(self.I)[:,axes_r]
        ls_r.I = np.array(list(set([tuple(i) for i in ls_r.I])))
                                        # Reinitialize active and old points
                                        # TODO: finish this... currently all
                                        # active.  Only needed for adaptivity
        ls_r.A, ls_r.O = set(range(len(ls_r.I))), set()
        return ls_r

    def add(self, k):
        """Add new active index k, array (dim,)."""
        self.I.append( np.array(k, dtype=np.int8) )
        self.A.add( len(self.I)-1 )
        return len(self.I)-1

    def active_to_old(self, i):
        """Convert active I-idx i to old index."""
        self.A.remove(i)
        self.O.add(i)    

    def get_active(self):
        return np.array(self.I)[list(self.A)]

    def get_all(self):
        return np.array(self.I)

    def max_level(self):
        return np.max(np.array(self.I))

    def plot(self, ax):
        """Plot LevelSet, each block is a level, with: black->active, 
        gray->old, white->not in set."""
        assert plt_loaded, 'Matplotlib not installed'
        assert self.dim == 2, 'Plotting only implemented in 2d'
        max_level = self.max_level()+1
        tlev = np.ones((max_level,)*self.dim, dtype=np.int8)*2
        for i in self.O:                # Old levels
            tlev[ tuple(self.I[i]-1) ] = 1   
        for i in self.A:                # Active levels
            tlev[ tuple(self.I[i]-1) ] = 0   
        ax.imshow(tlev, cmap=cm.gray, aspect='equal', origin='lower',
                  interpolation='nearest', 
                  extent=[.5,max_level+.5,.5,max_level+.5] )
        ax.get_xaxis().set_major_locator(plt.MaxNLocator(integer=True))
        ax.get_yaxis().set_major_locator(plt.MaxNLocator(integer=True))
        ax.set_xlabel('Level (dim 0)'); ax.set_ylabel('Level (dim 1)');


class SparseComponentRule:
    """
    A tensor-product rule forming a component of a sparse-grid rule.
    Exclusively for use in class SparseGrid.
    Members:
      dim      - dimension of rule - same as SparseGrid dimension
      delta    - False, standard rule U_i(f), True, \Delta = (U_i-U_{i-1})(f)
                 including interpolation (bool)
      quad     - nested quadrature object, e.g. QuadraturePatterson
      levelidx - quadrature rule level in each dim (list (dim,))
      nnodes   - number of nodes in each dimension (ndarray (dim,))
      i        - multi-indices of quadrature points, with reference to the finest
                 tensor-product grid, ndarray (N, dim) - N total points
      w        - weights of the multidimensional rule, corresponding to entries
                 of i, ndarray (N,)
    """
    def __init__(self, quad, levelidx, delta=False):
        i_1d = [quad.idx[l-1] for l in levelidx]
        w = quad.dw if delta else quad.w
        w_1d = [w[l-1] for l in levelidx]
        self.delta    = delta           # Original rule or difference rule
        self.levelidx = levelidx
        self.dim      = len(self.levelidx)
        self.quad     = quad            # 1d quadrature rule
        self.nnodes   = np.array([quad.nnode(l-1) for l in self.levelidx])
        self.i        = mylib.meshgrid_flatten(*i_1d)
        self.w        = mylib.meshprod_flatten(*w_1d)

    def nnodes_total(self): 
        """Total number of nodes in rule (scalar)"""
        assert self.w.shape[0] == np.prod(self.nnodes)
        assert len(self.w) == self.i.shape[0]
        return np.prod(self.nnodes)

    def global_index_to_x(self, midx):
        """Convert global-index to x-coordinates"""
        assert len(midx) == self.i.shape[1]
        return np.array(self.quad.global_index_to_x(np.array(midx)))

    def get_x(self):
        """Return all point coordinates as numpy array"""
        x = np.zeros((self.nnodes_total(), self.dim))
        for j, midx in enumerate(self.i):
            x[j,:] = self.global_index_to_x(midx)
        return x

    def extract_fvec(self, fval):
        """
        Given fval extract the values needed for this rule in a flat array
          fval - dictionary of function values indexed by tuples, function values
                 are allowed to be scalars or ndarrays (of any rank).
        Return:
          fvec - ndarray of function values, ordered according to self.i
        """
                                        # Detect datatype
        vex = fval.itervalues().next()  # Any entry from dict
        if type(vex) == np.ndarray:
            fvec = np.zeros((self.nnodes_total(),) + vex.shape)
        else:
            fvec = np.zeros(self.nnodes_total())
                                        # Fill array
        for j, midx in enumerate(self.i):
            loc = tuple(midx)
            assert fval.has_key(loc), 'Location not in fval: %s' % (loc,)
            fvec[j] = fval[loc]
        return fvec
        
    def lagrange(self, x):
        """
        Return all Lagrange polynomials for this rule at single location x,
        ordering according to self.i
        """
        assert len(x) == self.dim
        lag = self.quad.lagrange_delta if self.delta else self.quad.lagrange
        L_1d = [lag(x[j], self.levelidx[j]-1) for j in range(self.dim)]
        return mylib.meshprod_flatten(*L_1d)
     
    def integrate(self, fval):
        """
        Evaluate integral given data (dot-product of weights and values)
          fval - Either a) a dictionary of scalars/ndarrays indexed by tuples,
                 storage format of SparseGrid
                 Or b) an fvec - e.g. output of self.extract_fvec(fval)
        """
        fvec = self.extract_fvec(fval) if type(fval) == type({}) else fval
        return np.einsum('i...,i...', self.w, fvec)
        
    def interpolate(self, x, fval):
        """
        Polynomial interpolation based on this tensor-product rule
          x    - location at which to interpolate ndarray (dim,)
          fval - dictionary of function values indexed by tuples
        """
        return np.einsum('i...,i...', self.lagrange(x), self.extract_fvec(fval))

    def dimensionate(self, fvec):
        """
        Return fvec (itself the return value of self.extract_fvec()) viewed as a
        rank (self.dim) array - i.e. the most natural (unencoded) form of the
        data.  Nf is the size of the vector output in the case of non-scalar f.
        """
        if fvec.ndim == 1:              # Scalar-valued f
            return fvec.reshape(self.nnodes)
        else:                           # Tensor-valued f
            Nf = fvec.shape[1:]         # (shape of tensor-part)
            return fvec.reshape(np.hstack((self.nnodes, Nf)))
    
    def undimensionate(self, f):
        """Inverse of dimensionate - flatten f into fvec"""
        if f.ndim == self.dim:          # Scalar-valued f
            return f.reshape(self.nnodes_total())
        else:                           # Tensor-valued f
            Nf = f.shape[self.dim:]     # (shape of tensor-part)
            return f.reshape( (self.nnodes_total(),) + Nf )    

    def integrate_partial(self, fval, axes):
        """
        Integrate over a subset of the dimensions of the rule.  If axes lists all
        the axes this is equivalent to integrate().
          fval   - data evaluated at all support points of the component rule.  
                   If f is vector-valued the last dimension is the vector
          axes   - list of axes over which to integrate
        Return:
          sc_r   - SparseComponentRule object of reduced dimension
          fval_r - reduced rank version of fval, dictionary indexed by new
                   reduced-dim multi-indices
          axes_r - list of remaining axes in f_r (those not integrated yet)
        """
                                        # Full-rank representation of f
        f = self.dimensionate( self.extract_fvec(fval) )
        if len(axes) == 0: return f, range(self.dim)
                                        # Creative use of einsum()
        assert f.ndim < 26, ValueError('Implementation limit dim < 26')
        ijk  = np.array(list('abcdefghijklmnopqrstuvwxyz'))
                                        # Eg. cmd = 'abcde,b,c' for dim=3,
                                        # axis=[1,2], and matrix-valued f.
                                        # We avoid the ellipsis '...' because
                                        # of non-mature implementation.
        cmd  = ''.join(ijk[:f.ndim]) + ',' + ','.join(ijk[list(axes)])
        w    = self.quad.dw if self.delta else self.quad.w
        warg = [w[self.levelidx[ai]-1] for ai in axes]
        f_r = np.einsum(cmd, f, *warg)  # Partial integral of f
                                        # Complement of axes (remaining axes)
        axes_r = mylib.complement(axes, self.dim)
                                        # Remaining level indices, new comp rule
        levelidx_r = list(np.array(self.levelidx)[axes_r])
        sc_r = SparseComponentRule(self.quad, levelidx_r, delta=self.delta)
                                        # CHECK TODO: correct ordering?
        fvec_r = sc_r.undimensionate(f_r)
        fval_r = dict(zip([tuple(midx) for midx in sc_r.i], fvec_r))
        return sc_r, fval_r, axes_r


class SparseGrid:
    """
    Generalized sparse-grid with arbetrary (admissible) index sets.  Paper:
    Gerstner + Griebel, Dimension-Adaptive Tensor-Product Quadrature, Computing
    71(1), 2003.
    
    We use concept of "global-grid" for implementation - that is a full
    tensor-product grid in which all nodes of the sparse grid appear.  The size
    of the grid is determined by the maximum order rule in "quad".  Nodes of the
    sparse grid are refered to by their location in the full-grid, described by
    a multi-index (integer) array (dim,).  E.g. function values are stored in a
    dictionary with the corresponding multi-index as keys.  Therefore the
    efficiency of this implementation depends strongly on the performance of
    Python's native multi-index-keyed dictionaries.  Seems okay...

    Currently non-adaptive, but easy to extend.  Doesn't "flatten" weights etc.
    in the interest of preserving adaptive flexibility.

      dim        - dimension of the whole thing
      levelset   - LevelSet describing the sparse structure
      quad       - a nested quadrature object, e.g. QuadraturePatterson
      qcomponent - list of individual tensor-prod quadrature rules, of type
                   SparseComponentRule
    """
    def __init__(self, dim, quad, levelset=None, level=3, fill='simplex'):
        """
        Initialize either with a standard levelset (specify level and fill), or
        a custom levelset (specify levelset).  
        """
        self.dim      = dim
                                        # Init description of sparse pattern
        self.levelset = levelset if levelset else LevelSet(dim, level, fill=fill)   
        self.quad     = quad            # Init 1d rule
        self.init_tensorrules()         # Init component rules

    def init_tensorrules(self):
        """
        Init all tensor-product rule components of the sparse rule.  The
        levelset must be initialized first!  Store as list if
        SparseComponenetRule in qcomponent - ordering same as the ordering of I
        in LevelSet.
        """
                                        # Each index entry in levelset
                                        # represents an tensor-product
                                        # quadrature rule in dim
        I = self.levelset.get_all()     # dimensions.
        assert I.shape[1] == self.dim, 'Level set of wrong dimension'
        self.qcomponent = []
        for levelidx in I:
            self.qcomponent.append( SparseComponentRule(self.quad, levelidx,
                                                        delta=True) )
    def sample_fn(self, f, fargs=()):
        """
        Evaluate a function at all node locations in the sparse grid.  Return as
        fval, a dictionary of function values (possibly tensor) with keys tuples
        of node indices w.r.t. full grid.
        """
        fval = {}
        for q in self.qcomponent:
            for midx in q.i:
                loc = tuple(midx)
                if fval.has_key(loc): continue
                x = q.global_index_to_x(loc)
                fval[loc] = f(x, *fargs)
        return fval

    def integrate_partial(self, fval, axes):
        """
        Integrate over a subset of dimensions.  Return a new SparseGrid object
        defined over the remaining axes, together with new (integrated) sample
        values.
          fval - dictionary of measurements, e.g. sample_fn() output
          axes - dimensions over which to integrate, tuple of ints, non-
                 repeating.  If (0,...,dim-1) same as integrate()
        Return:
          sp    - a new SparseGrid class over the remaining axes
          fvalr - dictionary of measurements on the new instance sp
        """
        if len(axes) == 0:        return copy.deepcopy(self), fval
        if len(axes) == self.dim: return None, self.integrate(fval)
        dim_r = self.dim - len(axes)     # Remaining dimensions
                                         # Remaining levelset
        levelset_r = self.levelset.reduce_axes(axes)
                                         # New SparseGrid
        sp = SparseGrid(dim_r, self.quad, levelset=levelset_r)
        xval_r = sp.get_nodes()          # Nodes of reduced SparseGrid
                                         # Dict of new values initialized to 0.
        fval_r = dict.fromkeys(xval_r, 0.)
        for q in self.qcomponent:
                                         # Partial integration of component-rule
            qnew, fval_q, axes_r = q.integrate_partial(fval, axes)
                                         # Interpolation onto all new SparseGrid
                                         # nodes from the component rule
            for midx, x in xval_r.items():
                fval_r[midx] += qnew.interpolate(x, fval_q)
        return sp, fval_r

    def compute_sobol_indices(self, fval, cardinality=-1):
        """
        Compute Sobol indices based on Tang (2009), GLOBAL SENSITIVITY ANALYSIS 
        FOR STOCHASTIC COLLOCATION EXPANSION.  Approximates main-effect indices.
        Main effect indices (non-normalized) are returned in a dictionary D,
        with keys tuples of dimension indices.  E.g. D[(0,3,4)] is the
        main-effect index (cardinality 3) of the variables x_0, x_3, x_4.
        """
        cmax  = self.dim if cardinality==-1 else cardinality
        mu    = self.integrate(fval)    # Compute mean, variance
        fval2 = dict([(i, val**2) for i, val in fval.items()])
        var     = self.integrate(fval2) - mu**2
        D       = {(): mu**2}
        full    = range(self.dim)
        for c in range(1, cmax + 1):
                                        # Loop over all combinations of c indices
                                        # I.e. y = (0,), (1,), (dim-1,) then
                                        # y = (0,1), (0,2), (1,2), etc
            for y in itertools.combinations(full, c):
                z = mylib.complement(y, self.dim)
                                        # Integrate over dimensions z
                sp_r, fval_r = self.integrate_partial(fval, z)
                                        # Square result
                fval_r2 = dict([(i, val**2) for i, val in fval_r.items()])
                                        # Integrate over remaining dimensions
                D[y] = sp_r.integrate(fval_r2)
                for t in mylib.all_subsets(y, strict=True, null=True):
                    D[y] -= D[t]
        return D, mu, var

    def integrate(self, fval):
        return np.sum([q.integrate(fval) for q in self.qcomponent], axis=0)

    def interpolate(self, x, fval):
        return np.sum([q.interpolate(x, fval) for q in self.qcomponent], axis=0)

    #def sample_mc(self, fval, N):
    #    """Quasi-Monte-Carlo sampling of response: N - number of samples."""
    #    import sampling
    #    xs = sampling.sobol(N, self.dim)
    #    vex = fval.itervalues().next()  # Any entry from dict
    #    if type(vex) == np.ndarray: fs = np.zeros((N,) + vex.shape)
    #    else:                       fs = np.zeros(N)
    #    for i,x in enumerate(xs):
    #        if i % 1000 == 0:
    #            sys.stdout.write('%d '%i)
    #            sys.stdout.flush()
    #        fs[i] = self.interpolate(x, fval)
    #    return xs, fs
    
    def get_nodes(self):
        """
        Return spatial locations of nodes in dictionary (fval) format.
        Use: np.array(get_nodes().values()) to get flat-array format.
        """
        return self.sample_fn(lambda x: x)

    def get_weights(self):
        """Return weights of the rule in dictionary (fval) format."""
        wval = {}
        for q in self.qcomponent:
            for midx, w1 in zip(q.i, q.w):
                if tuple(midx) in wval: wval[tuple(midx)] += w1
                else:                   wval[tuple(midx)]  = w1
        return wval

    def n_nodes(self):
        return len(self.get_nodes())

    def plot(self, outfile='a.png'):
        """Plot the levelset and sparse grid."""
        assert plt_loaded, 'Matplotlib not installed'
        assert self.dim==2 or self.dim==3, 'Plotting only implemented in 2d/3d'
        if self.dim == 2:
            fig = plt.figure(figsize=(10,5))
            ax = fig.add_subplot(121); self.levelset.plot(ax)
            ax = fig.add_subplot(122)
            x = np.array(self.get_nodes().values())
            ax.plot(x[:,0], x[:,1], 'xk')
            fig.savefig(outfile)
        elif self.dim == 3:
            x = np.array(self.get_nodes().values())
            from mpl_toolkits.mplot3d import Axes3D
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(x[:,0], x[:,1], x[:,2])
            fig.show()


def integrate_unitcube(dim, level, f, fargs=()):
    """
    Wrapper to do integration on [0,1]^dim with class SparseGrid.  Example of
    basic use of SparseGrid class.
    """
    sp = SparseGrid(dim, QuadraturePatterson(),
                    level=level, fill='simplex')
    fval = sp.sample_fn(lambda x: f(x, *fargs))
    return sp.integrate(fval)

def sparse_nodes(dim, level, fill='simplex'):
    """
    Extract nodes and weights of a sparse rule of given dimension (dim), level
    and fill type.
    Return:
      midx - list of node multi-indices, needed for constructing dictionary
             needed for calling SparseGrid.integrate()
      x    - list of node locations
      w    - list of weights
    """
    sp = SparseGrid(dim, QuadraturePatterson(),
                    level=level, fill=fill)
    xval, wval = sp.get_nodes(), sp.get_weights()
    midx, x = xval.keys(), xval.values()
                                        # Guarantee weights are in the same order
    w = np.zeros(len(midx))             # as the nodes and multi-indices.
    for i, mi in enumerate(midx): w[i] = wval[mi]
    return sp, midx, x, w


#===============================================================================
# Unit tests
#===============================================================================
def unittest_integration(fill='simplex', quadrature=QuadraturePatterson()):
    """Verify integration implementation against Genz functions."""
    def const(x): return 1.
    import test_genz as genz

    for dim in [2,3,4,5,8,10,12,16,20]:
        def genz_mod(x):         return genz.f1(x)
        def genz_exact_mod(dim): return genz.f1_exact(dim)

        genz.c = np.array([0.5] * dim); genz.c = genz.c / sum(genz.c) * 9.
        genz.w = np.array([0.5] * dim)
        exact = genz_exact_mod(dim)
        print '%s Dimension = %d %s' % ('-'*20, dim, '-'*20)
        print '%10s %10s %16s %12s' % ('Level','#nodes','Integral','Rel error')
        for l in range(1, quadrature.nlevel()+1):
            sp = SparseGrid(dim, quadrature, level=l, fill=fill)
            if dim==2: sp.plot('sparse_simplex_%d.png' % l)
            fval = sp.sample_fn(genz_mod)
            approx = sp.integrate(fval)
            print '%10d %10d %16.10g %12.1e' % (l, len(fval), approx,
                                                abs(approx-exact)/exact)
        print '%10s %10s %20.10g' % ('Exact','',exact)

def unittest_interpolation(fill='simplex', quadrature=QuadraturePatterson()):
    """Verify interpolation implementation against some function f().
    Plot key - black surface: original, red: reconstruction."""
    dim, level = 2, 5
    sp = SparseGrid(dim, quadrature, level=level, fill=fill)
    def f(x): return np.cos(2*np.pi*x[0])*np.cos(2*np.pi*x[1]+2)
    fval = sp.sample_fn(f)
    sx   = np.array(sp.get_nodes().values())
                                    # Plotting interpolation
    M    = 41
    x1   = np.linspace(0,1,M)
    xy   = mylib.meshgrid_flatten(x1, x1)
    X, Y = xy[:,0].reshape(M,M), xy[:,1].reshape(M,M)
                                    # 2d sampling
    F, Fa   = np.zeros(xy.shape[0]), np.zeros(xy.shape[0])
    for i,x in enumerate(xy):
        F[i]  = f(x)
        Fa[i] = sp.interpolate(x, fval)
                                    # Plotting
    if plt_loaded:
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure(figsize=(12,8))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_wireframe(X, Y, F.reshape(M,M), colors='k')
        ax.plot_wireframe(X, Y, Fa.reshape(M,M), colors='r')
        ax.scatter(sx[:,0], sx[:,1], 0, c='k') 
        fig.savefig('unittest_interpolation.pdf')

def unittest_sobol(fill='simplex', quadrature=QuadraturePatterson()):
    """Verify Sobol index implementation against simple polynomial."""
    dim = 2
    def integrand(x):                   # Rosenbrock function
        x, y = tuple(x)
        return 100.0*(y-x*x)*(y-x*x) + (1.-x)*(1.-x)
                                        # Sparse-grid computation
    sp = SparseGrid(dim, quadrature, level=4, fill=fill)
    fval = sp.sample_fn(integrand)
    D, mu, var = sp.compute_sobol_indices(fval, cardinality=-1)
    sys.stdout.write('Sobol indices (non-normalized main-effect indices)\n')
    sys.stdout.write('  Sparse grid results\n')
    sys.stdout.write('      Mean  = %10.4g\n' % mu)
    sys.stdout.write('      Var   = %10.4g\n' % var)
    for i in range(dim):
        sys.stdout.write('      D_%1d   = %10.4g\n' % (i, D[(i,)]))
    for i in range(1,dim):
        for j in range(i):
            sys.stdout.write('      D_%1d,%1d = %10.4g\n' % (j,i,D[(j,i)]))
                                        # Reference calculations -
                                        # requires module sobol_mc
    import sobol_mc
    sys.stdout.write('  Monte-Carlo results\n')
    for i in range(dim):
        M, mu, var, D, D_tot = sobol_mc.sobol_variance_mc(integrand, dim,
                                                [i], N=100000, monitor=False)
        if i == 0:
            sys.stdout.write('      Mean  = %10.4g\n' % mu)
            sys.stdout.write('      Var   = %10.4g\n' % var)
        sys.stdout.write('      D_%1d   = %10.4g\n' % (i, D))
            
def unittest_nodecount(fill='simplex', quadrature=QuadraturePatterson()):
    """Number of nodes required for various dimensions and levels - pretty-print.
    This will depend on the fill-type and quadrature rule.  Handy for choosing a
    level."""
    dims = [2,3,4,5,7,8,10,12,16,20]
    levels = range(1, quadrature.nlevel()+1)
    sys.stdout.write('Table: Number of support-points in sparse grid\n')
    sys.stdout.write('Fill type: %s\n' % fill)
    sys.stdout.write(' '*40+'level\n             ')
    for level in levels: sys.stdout.write('%10d' % level)
    sys.stdout.write('\n            '+'-'*6*11)
    for dim in dims:
        sys.stdout.write('\n %s %6d |' % ('dim' if dim==8 else '   ', dim))
        for level in levels:
            sp = SparseGrid(dim, quadrature, level=level, fill=fill)
            sys.stdout.write('%10d' % sp.n_nodes()); sys.stdout.flush()
    sys.stdout.write('\n')


if __name__ == '__main__':

    if True:                            ### Test: count support-points
        unittest_nodecount(fill='simplex', quadrature=QuadratureCC())

    if False:                           ### Test: integration
        unittest_integration(fill='simplex', quadrature=QuadraturePatterson())
    
    if False:                           ### Test: Sobol indices
        unittest_sobol(fill='simplex', quadrature=QuadratureCC())

    if False:                           ### Test: interpolation
        unittest_interpolation(fill='simplex', quadrature=QuadratureCC())

    if False:                           ### Example: typical use with cheap fn
        dim, level = 4, 3
        def fcheap(x): return np.sin(x[0])* x[1] + x[2]**2 * x[3]
            
        sp = SparseGrid(dim, QuadratureCC(), level=level, fill='simplex')
        fval = sp.sample_fn(fcheap)
        print 'Integral via SparseGrid: %10.4g' % sp.integrate(fval)
        
    if False:                           ### Example: typical use with expensive
                                        ### fn requiring external evaluation
 
        dim, level = 4, 3               # Get multi-indices, nodes and weights
        sp, midx, x, w = sparse_nodes(dim, level, fill='simplex')
                                        # Sample function at nodes - this part
                                        # can be split off including backup,
                                        # error handling, parallization etc.
        def fexpensive(x): return np.sin(x[0])* x[1] + x[2]**2 * x[3]
        fi = np.zeros(len(w))
        for i, xi in enumerate(x): fi[i] = fexpensive(xi)
                                        # Integral via dot product...
        print '    Integral via dot-product: %10.4g' % np.dot(fi, w)
                                        # OR reconstruct fval, and apply
                                        # SparseGrid.  This is the only way to
                                        # compute Sobol indices etc.
        fval = dict(zip(midx, fi))
        print 'Integral via SparseGrid call: %10.4g' % sp.integrate(fval)
