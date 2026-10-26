"""Whole-interval positivity check for an unnormalised Chebyshev background."""
import numpy as np
from numpy.polynomial import Chebyshev

def chebyshev_minimum(coefficients):
    """Minimum of 1 + sum(a_i T_(i+1)) on [-1, 1], including all extrema."""
    c=np.asarray([1.,*coefficients],dtype=float)
    if not np.isfinite(c).all(): raise ValueError('Nonfinite Chebyshev coefficient')
    p=Chebyshev(c)
    points=[-1.,1.]
    for z in p.deriv().roots():
        if abs(z.imag)<1e-10 and -1.<z.real<1.: points.append(float(z.real))
    return float(np.min(p(np.asarray(points,dtype=float))))
