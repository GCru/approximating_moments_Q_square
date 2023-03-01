from math import atan, exp, sin
from scipy.integrate import quad
import numpy

from momentchi2 import hbe, lpb4,sw,wf

def theta(u, lambd, h, x, delta2):
    """
    See Imhof (1961), p.423
    """
    m = len(lambd)
    sum_val = 0.0
    for i in range(m):
        sum_val += h[i] * atan(lambd[i] * u) + delta2[i] * lambd[i] * u / (1.0 + (lambd[i] * u)**2)
    sum_val = 0.5 * sum_val - 0.5 * x * u
    return sum_val

def rho(u, lambd, h, delta2):
    """
    See Imhof (1961), p.423
    """
    m = len(lambd)
    prod = 1.0
    for i in range(m):
        prod *= (1.0 + (lambd[i] * u)**2)**(0.25 * h[i]) * exp(0.5 * delta2[i] * (lambd[i] * u)**2 / (1.0 + (lambd[i] * u)**2))
    return prod

def imhoffunc(u, lambd, h, x, delta2):
    """
    This is the function under the integral sign in equation (3.2), Imhof (1961), p.422
    """
    res = sin(theta(u, lambd, h, x, delta2)) / (u * rho(u, lambd, h, delta2))
    return res


def imhoff(x, lambd, h, delta2,  epsabs=1.49e-05, epsrel=1.49e-05, limit=10000):
    """
    Implements the Imhof (1961) algorithm
    """
   
    Qx = 0.5 + quad(imhoffunc, 0, numpy.inf, args=(lambd, h, x, delta2,), epsabs=epsabs, epsrel=epsrel, limit=limit)[0]/numpy.pi

    return Qx

if __name__ == '__main__':
    
    coeff = [0.1]*10
 
    x=1.0
    
    print(imhoff(x,coeff,[1]*len(coeff),[0]*len(coeff)))

    print(1-lpb4(coeff=coeff, x=x, p=4))
