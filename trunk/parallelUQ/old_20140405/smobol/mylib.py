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
"""
Module of useful generic functions.

My own python routines for relatively generic tasks.  Occasionally I
discover that one of these already exists in the standard libraries -
in which case the version here is depreciated.  Already depreciated are:

* Remove infs and nans (for plotting) - use numpy.nan_to_num()
* Norm - use numpy.linalg.norm
* Linspace like in MatLAB - use numpy.linspace.

Numpy and Scipy are missing some strange things however - like forward_sub()
for solving a lower-triangular linear system.  Unless otherwise stated,
tensor arguments are expected to be numpy arrays.
"""
import numpy as np
import os, sys, re, shutil, copy, time
import cPickle, itertools

import mylib_meshgrid
#---------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------
def pprint_array(a):
    m,n = a.shape
    for i in range(m):
        sys.stdout.write('[')
        for j in range(n):
            sys.stdout.write('%12.6g' % a[i,j])
        sys.stdout.write(']\n')

def timestamp():
    return time.strftime('%Y-%m-%d_%H%M%S')

#---------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------
def fnname():
    """Return name of the calling function."""
    return sys._getframe(1).f_code.co_name + '()'

def fncountcall(function, args):
    """
    Count the number of calls of a function.  Call once at the start of a
    routine, and henceforth use the returned function instead of the original.
    The other return value will contain the number of function calls updated.
    (Trick from scipy.optimize.)
    """
    ncalls = [0]
    def fn_wrapper(x):
        ncalls[0] += 1
        return function(x, *args)
    return fn_wrapper, ncalls

def flatten(x):
    """
    Returns a single, flat list which contains all elements retrieved
    from the sequence and all recursively contained sub-sequences (iterables).
    For example: ::

      >>> flatten([[[1,2,3], (42,None)], [4,5], [6], 7, MyVector(8,9,10)])
      [1, 2, 3, 42, None, 4, 5, 6, 7, 8, 9, 10]
    """
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result

#---------------------------------------------------------------------------
# String manipulations.
#---------------------------------------------------------------------------
def str2list(str):
    """Convert string to a list of strings."""
    return str.split('\n')

def list2str_str(list):
    """Convert list to a single string.  Inverse of str2list."""
    return ''.join([str(i)+'\n' for i in list])

def list2str(list, format='%13.5e'):
    """More general version of list2str_str."""
    out = ''
    for i in tolist(list):
        out = out + format % (i,)
    return out

def read_float_array(filename, skipline=0):
    """
    Read an array of floats into an a list-of-lists.  Skip skipline
    lines from the start of the file.  Don't convert to a numpy array
    - to allow for varying numbers of entries per line.
    """
    a = open('data_raw2.txt', 'r').readlines()[skipline:]
    return [[float(j) for j in i.split()] for i in a]

#---------------------------------------------------------------------------
# Regular expressions on files - Unix-ey way of thinking led to these...
#---------------------------------------------------------------------------
def re_infile(filename, restr, group):
    """
    Get all matches of a regular expression in a file.  Return only the group
    specified in each match, in the form of a list.
    """
    if not os.path.isfile(filename):
        print >> sys.stderr, 're_infile(): File not found.'
        return None
    value = []
    fh = open(filename, "r")
    reobj = re.compile(restr)
    for line in fh:
        m = reobj.match(line)  #search
        if m: value.append(m.group(group))
    fh.close()
    return value


def re_infile_last(filename, restr, group):
    """
    Get the *last* match of a regular expression in a file.  Return only the
    group specified.
    """
    tmp = re_infile(filename, restr, group)
    if len(tmp) > 0: return tmp[-1]
    else:            return None


def re_infile_replace(filename, oldstr, newstr):
    """
    Replace every instance of a re in a file with another string.  Return
    number of replacements performed.  Uses temporary file :(
    """
    nreplace = 0
    fh = open(filename, "r")
    fhw = open(filename + ".pytmp", "w")
    
    reobj = re.compile(oldstr)
    for line in fh:
         (repstr, n) = reobj.subn(newstr, line)
         print >> fhw, repstr,
         nreplace = nreplace + n
    shutil.move(filename + ".pytmp", filename)
    return nreplace
        

#---------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------
def tolist(a):
    """Given a list or a scalar, return a list."""
    if getattr(a, '__iter__', False): return a
    else:                             return [a]

def toarray(a):
    """Given an array, list, or scalar, return a np.array."""
    if getattr(a, '__iter__', False): return np.array(a)    # List case
    else:                             return np.array([a])  # Scalar case

def tomat(a):
    """
    Arg is an numpy array.  If it's a vector, return corresponding n*1
    matrix.  Otherwise return arg unchanged.
    """
    if a.ndim == 1:   return np.array([a]).T
    elif a.ndim == 2: return a
    else:             raise ValueError("Arg must be a numpy array of ndim=1,2.")

def mylen(a):
    """
    If a is a scalar (or object with no length) return 1, otherwise
    len(a).  Slight generalization of built-in len().
    """
    if '__len__' in dir(a): return len(a)
    else:                   return 1

def floateq(a, b, eps=1.e-12):
    """Are two floats equal up to a difference eps?"""
    return abs(a - b) < eps

def project_normal(x, c1):
    """Project vector x onto plane normal to vector c1."""
    c = c1 / np.linalg.norm(c1, ord=2)
    return x - c * np.dot(x, c)

#---------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------
def complement(y, n):
    """Complement of y in [0,...,n-1]. E.g. y = (0,2), n=4, return [1,3]."""
    return sorted(list(set(range(n)) - set(y)))            

def all_subsets(A, strict=True, null=True):
    """
    Return all subsets of A (list/tuple/etc).  If strict is False, result
    includes A, if null is True, includes ().
    """
    n = len(A)
    min = 0 if null else 1
    max = n if strict else n+1
    out = []
    for m in range(min, max):
        out.append( itertools.combinations(A, m) )
    return itertools.chain(*out)



#---------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------
def interpolate_zeros(x, f):
    """
    For two float arrays of data, representing a 2d plot (x, f) in
    order - determine all the locations of x at which f crosses zero
    (using linear interpolation).  Return as a list.
    """
    assert len(x) == len(f)
    xcross = []
    for i in range(len(f) - 1):
        x0, x1 = x[i], x[i+1]
        f0, f1 = f[i], f[i+1]
        if (f0 >= 0. and f1 < 0.) or (f0 <= 0. and f1 > 0.):
            xcross.append(x0 - (x1-x0) * f0/(f1-f0))
    return xcross


#---------------------------------------------------------------------------
# Solve lower/upper triangular systems apparently nowhere in numpy/scipy.
#---------------------------------------------------------------------------
def forward_sub_single(A, b):
    """Solve Ax = b for A lower triangular for single RHS b."""
    if A.shape[0] != b.shape[0]:
        raise ValueError ("A and b must be compatible.")
    x = np.zeros(A.shape[1])
    for j, row in enumerate(A):
        pivot = row[j]
        if floateq(pivot, 0):
            raise ValueError ("Matrix has a zero diagonal element.")
        x[j] = (b[j] - np.dot(row, x)) / pivot
    return x

def backward_sub_single(A, b):
    """Solve Ax = b for A upper triangular for single RHS b."""
    if A.shape[0] != b.shape[0]:
        raise ValueError ("A and b must be compatible.")
    x = np.zeros(A.shape[1])
    rev = range(A.shape[0]); rev.reverse()
    for j in rev:
        row = A[j]
        pivot = row[j]
        if floateq(pivot, 0):
            raise ValueError ("Matrix has a zero diagonal element.")
        x[j] = (b[j] - np.dot(row, x)) / pivot
    return x

def tri_sub_generic(A, b, sub_single_fn):
    """Helper fn for forward_sub() and backward_sub(), enabling multiple
    RHSs in both cases."""
    if b.ndim == 1:
        return sub_single_fn(A, b)
    else:
        outp = np.zeros(b.shape)
        for i in range(b.shape[1]):
            outp[:,i] = sub_single_fn(A, b[:,i])
        return outp

def forward_sub(A, b):
    """Solve Ax = b for A lower triangular, possible multiple RHSs."""
    return tri_sub_generic(A, b, forward_sub_single)

def backward_sub(A, b):
    """Solve Ax = b for A upper triangular, possible multiple RHSs."""
    return tri_sub_generic(A, b, backward_sub_single)


def complex_derivative(fn, x, eps=1.e-12):
    """
    Compute the 1st-derivatives of a real fn of a single array
    argument (x), returning a single array (f) using complex
    differences, about x.  If the function has len(x)=N arguments,
    this requires N *complex* calls of the function, costing typically
    2N real calls.  Insensitive to epsilon.  Return an M x N Jacobian
    matrix.  TODO: case handling for funcs with scalar input/output?
    """
    assert x.ndim == 1
    N = x.shape[0]
    # Calculate the 1st-difference to discover the shape of the output.
    df0 = fn(x + unitvec(0, N) * eps * 1.j).imag / eps
    assert df0.ndim == 1
    df = np.zeros((df0.shape[0], N))
    df[:,0] = df0
    # Calculate remaining differences.
    for i in range(1, N):
        df[:,i] = fn(x + unitvec(i, N) * eps * 1.j).imag / eps
    return df

def feuler(f, t, y, dt):
    """
    Forward-Euler simplest possible ODE solver - 1 step.

    :param f:  function defining ODE dy/dt = f(t,y), two arguments
    :param t:  time
    :param y:  solution at t
    :param dt: step size
    :return: approximation of y at t+dt.
    """
    return y + dt * f(t, y)

def rk4(f, t, y, dt):
    """
    Runge-Kutta explicit 4-stage scheme - 1 step.

    :param f:  function defining ODE dy/dt = f(t,y), two arguments
    :param t:  time
    :param y:  solution at t
    :param dt: step size
    :return: approximation of y at t+dt.
    """
    k1 = f(t, y)
    k2 = f(t + dt/2, y + dt*k1/2)
    k3 = f(t + dt/2, y + dt*k2/2)
    k4 = f(t + dt, y + dt*k3) 
    return y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)

#---------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------
def cond(A):
    """Find condition number of a matrix by SVD (brute force)."""
    s = np.linalg.svd(A)[1]
    return s[0] / s[-1]
    
def shift_array(l, offset):
    """Shift an array l by offset."""
    offset %= len(l)
    return np.concatenate((l[-offset:], l[:-offset]))

def gridarray(a, b):
    """
    Given two arrays create an array of all possible pairs, a 2d grid.
    E.g. a = [1, 2], b = [2, 4, 5], gridarray(a,b) = [[1,2], [1,4],
    [1,5], [2,2], [2,4], [2,5]].  May be used repeatedly for increasing
    dimensionality.
    DEPRECIATED: Use A, B = np.meshgrid(a, b).
                 Note that meshgrid works with arbetrary dimension too.
    """
    if a == None: return b   # Trivial cases
    if b == None: return a   
    adim, bdim = 1, 1
    if a.ndim > 1: adim = a.shape[1]
    if b.ndim > 1: bdim = b.shape[1]
    ab = np.zeros((a.shape[0] * b.shape[0], adim + bdim),dtype=a.dtype)
    count = 0
    for aa in a:
        for bb in b:
            ab[count,0:adim] = aa
            ab[count,adim:]  = bb
            count = count + 1
    return ab

def meshgrid_flatten(*X):
    """
    Functionally same as numpy.meshgrid() with different output
    format.  Function np.meshgrid() takes n 1d ndarrays of size
    N_1,...,N_n, and returns X_1,...,X_n n-dimensional arrays of shape
    (N_1, N_2,... N_n).  This returns instead a 2d array of shape
    (N_1*...*N_n, n).
    """
    if len(X) == 1:                     # Because np.meshgrid() can't handle
        return np.array([X[0]]).T       # less than 2 arguments
    return np.vstack( map(lambda x: x.flatten(),
                          mylib_meshgrid.meshgrid(*X, indexing='ij')) ).T

def meshprod_flatten(*X):
    """
    Given n arrays create an array of all combinations of products.
    E.g.: a = [1, 2]; b = [2, 4, 5]
          meshprod_flatten(a,b) = [2,4,5,4,8,10]
    in the same order as the result of meshgrid_flatten().  This is
    useful for e.g. constructing tensor products of quadrature
    weights.
    """
    return np.prod(meshgrid_flatten(*X), axis=1)

def gridprod(a, b):
    """RPD: Depreciated, use meshprod_flatten()."""
    return (np.array([a]).T * np.array([b])).flatten()

def reverse(list):
    """Reverse a list out-of-place."""
    l = copy.deepcopy(list)
    l.reverse()
    return l

def gridlist(a, b):
    """
    As for gridarray(), but for Python lists - only 1d input lists allowed.
    """
    ab = []
    for aa in a:
        for bb in b:
            ab.append([aa, bb])
    return ab

#---------------------------------------------------------------------------
# Dictionary operations.
#---------------------------------------------------------------------------
class Bunch(object):
    """
    If you have a dictionary d and want to access (read) its values
    with the syntax x.foo instead of the clumsier d['foo'], just do
    x = Bunch(d).  TODO: what about write - that works too, right?
    """
    def __init__(self, adict):
        self.__dict__.update(adict)


def find_key(dic, val):
    """Return the keys of dictionary dic with the value val."""
    return [k for k, v in dic.iteritems() if v == val]

def dict_subset(dic, keys):
    """Extract subset of keys from dictionary."""
    return dict((k, dic[k]) for k in keys if k in dic)

#---------------------------------------------------------------------------
# Sorting routines.
#---------------------------------------------------------------------------
def sort(list):
    """Sort and return."""
    a = copy.deepcopy(list)
    a.sort()
    return a

def sort_index(list):
    """
    Return the *indices* of the elements in an array in order from smallest
    element to largest, e.g sort_index([3, -1, 1]) gives [1, 2, 0].
    """
    x = zip(list, range(0, len(list)))
    x = sorted(x, key=lambda x: x[0])
    return map(lambda x: x[1], x)

def minindex(A):
    """Return the index of the minimum entry in a list.  If there are 
    multiple minima return one."""
    return min((a, i) for i, a in enumerate(A))[1]

def minindices(A):
    """
    Return a list of indices of the (possibly multiple) minimum
    elements in list of floats A.
    """
    I = sort_index(A)
    Imin = [I[0]]
    for i in range(np.size(I) - 1):
        if not floateq(A[I[i]], A[I[i+1]]): break
        Imin.append(I[i+1])
    return Imin

def reorder_index(list, imap):
    """
    Given an index mapping - e.g. from sort_index() - reorder a list
    according to this mapping.  This is what numpy arrays do
    automatically with the syntax a[imap].
    TODO: implement this without deepcopy.
    """
    tmp = copy.deepcopy(list)
    for i, j in enumerate(imap): tmp[i] = list[j]
    return tmp

def sort_list2(a, b):
    """
    Sort list b according to cmp applied to a.  Return both a and b.
    E.g.  if a = [2, 0, 1] and b = ['a', 'b', 'c'] return [0, 1, 2],
    ['b', 'c', 'a'].
    """
    z = zip(*sorted(zip(a, b), key=lambda x: x[0]))
    return list(z[0]), list(z[1])

def listmatch(a, b):
    """
    True if a contains same float entries as b where b is not None.
    E.g. listmatch([1.,2.,3.], [1., None, 3.]) returns True.
    """
    if not len(a) == len(b):
        raise ValueError("Arguments must be same length, a=%s, b=%s."
                         % (str(a), str(b)))
    for ai, bi in zip(a,b):
        if bi is not None:
            if not floateq(ai, bi):
                return False
    return True

def closest(x, a):
    """
    Find the index of the element of the float list a closest to x.
    """
    min_i, min_dist = -1, 1.e99
    for i, y in enumerate(a):
        if fabs(x - y) < min_dist:
            min_i, min_dist = i, fabs(x - y)
    return min_i

def interval_idx(x, xi):
    """
    Given a grid of points xi (sorted smallest to largest) find in
    which interval the point x sits.  Returns i, where x \in [ xi[i],xi[i+1] ].
    Raise ValueError if x is not inside the grid.
    """
    assert x >= xi[0] and x <= xi[-1], ValueError('x=%g not in [%g,%g]' % (x,xi[0],xi[-1]))
    for i in range(len(xi)-1):
        if xi[i] <= x and x <= xi[i+1]:
            return i

#---------------------------------------------------------------------------
# Matlab functions.
#---------------------------------------------------------------------------
def linspace(min, max, n):
    """
    Return n uniform samples in interval [min, max] including
    endpoints.  Equivalent to Matlab's linspace().  V2: discovered
    numpy has the function, DEPRECIATED.
    """
    return np.linspace(min, max, n)

def unitvec(m, n):
    """
    Return the m-th Cartesian basis vector in n-dimensions.
    """
    tmp = np.zeros(n)
    tmp[m] = 1.
    return tmp

def symmetric_banded_matrix(x):
    """
    Create a symmetric banded matrix of dimension NxN from a vector x
    of dimension N.  The diagonal values are x[0], the upper and lower
    off-diagonal x[1], etc.
    """
    N = x.shape[0]
    m = np.zeros((N, N))
    for i in range(N):
        for j in range(i+1):
            m[i,j] = m[j,i] = x[i-j]
    return m

def fzero_repeated_bisection(fn, x0=[0.,1.], max_iter=100, x_tol=1.e-8):
    """
    Solve fn(x) = 0 using repeated bisection - 1d only.  Initial guess and
    result are both intervals.

    :param fn:       Function to solve.
    :param x0:       Initial interval [a, b], must contain a root.
    :param max_iter: Max steps, final interval is 1/2**max_iter times smaller
                     than initial interval.
    :param x_tol:    Terminate when interval width < tol.

    :return:         Final interval x.
    """
    x = copy.copy(x0)
    ftmp = [fn(x[0]), fn(x[1])]  # Storage: never calculate fn twice at same location.
    if ftmp[0]*ftmp[1] > 0:
        raise ValueError("Function value doesn't change sign on initial interval: "
                         "[%e, %e], f = [%e, %e]." % (x[0],x[1],ftmp[0],ftmp[1]))
    for i in range(max_iter):
        midpnt = (x[1] + x[0]) / 2.
        fmid   = fn(midpnt)
        if ftmp[0]*fmid < 0:
            x[1], ftmp[1] = midpnt, fmid
        else:
            x[0], ftmp[0] = midpnt, fmid
        if abs(x[1] - x[0]) < x_tol: break
    return x



def linecut(x, y, ycut):
    """
    Given a sequence of line segments defined by arrays x, y, and a
    specified y value ycut, return a list of x values where the line
    segments cross ycut, return empty list if they don't cut.  

    NB: Often useful when e.g. x is an iteration number and y a
    residual, xcut gives then at what iterations a convergence
    criteria was achieved.
    """
    xcut = []
    if x.size < 2: return []
    for i in range(x.size-1):
        xx, yy = [x[i], x[i+1]], [y[i], y[i+1]]
        if yy[1] < yy[0]: 
            xx.reverse()
            yy.reverse()
        if ycut >= yy[0] and ycut < yy[1]:
            yfrac = (ycut - yy[0]) / (yy[1] - yy[0])
            xcut.append(xx[0] + yfrac * (xx[1] - xx[0]))
    return xcut


def linrecon_cartesian(coord, samp, grad, grad_eps=1.e-3):
    """
    Given a single coordinate, value and gradient, return a list of 
    coordinates and values obtained by linear reconstruction in all 
    Cartesian directions (not including the original point).
    """
    ndim = coord.shape[0]        
    naug = 2*ndim
    augcoord = np.zeros((naug, ndim))
    augsamp  = np.zeros(naug)
    for j in range(ndim):
        dx = np.zeros(ndim); dx[j] = 1.
        augcoord[2*j]   = coord + grad_eps * dx
        augcoord[2*j+1] = coord - grad_eps * dx
        augsamp[2*j]    = samp  + grad_eps * grad[j]
        augsamp[2*j+1]  = samp  - grad_eps * grad[j]
    return augcoord, augsamp


#---------------------------------------------------------------------------
# Basic numerics
#---------------------------------------------------------------------------
def cubic_spline(xx, xi, fi):
    """
    Cubic spline in 1d.
      xx - prediction locations
      xi - node locations, ordered
      fi - node values
    Return - prediction values.  Will fail with ValueError if one of the xx
    is outside the interval [xi[0], xi[-1]].
    """
    def S(xp, xj, a, b, c, d):
        return a + b * (xp-xj) + c * (xp-xj)**2 + d * (xp-xj)**3
    h  = xi[1:] - xi[:-1]                         # (n,)
    hh = 2*(h[1:] + h[:-1])                       # (n-1,)
    n  = len(h)
    A = np.zeros((n+1,n+1))                       # Matrix
    for i in range(1, n):
        A[i,i], A[i,i-1], A[i,i+1] = hh[i-1], h[i-1], h[i]
    A[0,0], A[n,n] = 1, 1
    rhs = np.zeros(n+1)                           # RHS
    for i in range(1, n):
        rhs[i] = 3.*(fi[i+1] - fi[i])/h[i] - 3.*(fi[i] - fi[i-1])/h[i-1]
    c = np.dot(np.linalg.inv(A), rhs)             # Solve system - (n+1,)
                                                  # Get remaining coeffs
    a = fi[:]                                                  # (n+1,)
    b = (a[1:] - a[:-1]) / h[:] - h[:]/3. * (2*c[:-1] + c[1:]) # (n,)
    d = (c[1:] - c[:-1]) / (3. * h[:])                         # (n,)
    out = np.zeros(xx.shape[0])                   # Prediction
    for i,x in enumerate(xx):
        idx = interval_idx(x, xi)
        out[i] = S(x, xi[idx], a[idx], b[idx], c[idx], d[idx])
    return out

def linear_spline(xx, xi, fi):
    """
    Linear spline in 1d - i.e. linear interpolation of nodes.
      xi - node locations, ordered
      fi - node values
      xx - prediction locations
    Return - prediction values
    """
    def S(xp, xj, a, b): return a + b * (xp-xj)
    h  = xi[1:] - xi[:-1]                         # (n,)
    n  = len(h)
    a = fi[:]                                     # (n+1,)
    b = (fi[1:]-fi[:-1]) / h[:]
    out = np.zeros(xx.shape[0])
    for i,x in enumerate(xx):
        idx = interval_idx(x, xi)
        out[i] = S(x, xi[idx], a[idx], b[idx])
    return out



#-------------------------------------------------------------------------
# Unit test.  Because there is so much independent stuff in this
# module, this is only a test of some more complex routines.
#-------------------------------------------------------------------------
if __name__ == '__main__':

    # Triangular matrix solution A x = b.
    dim = 10
    A = np.array(range(1, dim**2+1)).reshape((dim,dim))
    Aupper, Alower = copy.copy(A), copy.copy(A)
    b = np.array(range(dim))
    
    for i in range(dim):
        for j in range(i):
            Aupper[i][j] = 0.
            Alower[j][i] = 0.

    xlower = forward_sub(Alower, b)
    xupper = backward_sub(Aupper, b)
    
    blower = np.tensordot(Alower, xlower, axes=(1,0))
    bupper = np.tensordot(Aupper, xupper, axes=(1,0))

    for i in range(dim):
        print '%16.6f %16.6f %16.6f %16.6f' % \
            (xlower[i], xupper[i], blower[i], bupper[i])
