import numpy as np

def spherical_harmonic(m, l, theta, phi):
    """
    Compute the spherical harmonic of degree l and order m
    at the angles theta and phi using recursion relations.
    """
    if m == 0:
        return np.sqrt((2*l+1)/(4*np.pi)) * legendre(l, 0, np.cos(theta)) * np.exp(1j*m*phi)
    elif m < 0:
        return (-1)**m * np.sqrt((2*l+1)/(4*np.pi)) * legendre(l, -m, np.cos(theta)) * np.exp(1j*m*phi)
    else:
        return np.sqrt((2*l+1)/(4*np.pi)) * legendre(l, m, np.cos(theta)) * np.exp(1j*m*phi)

def legendre(l, m, x):
    """
    Compute the associated Legendre polynomial of degree l and order m
    at the point x using recursion relations.
    """
    if m == 0:
        return np.polynomial.legendre.Legendre.basis(deg=l)(x)
    elif m < 0:
        return (-1)**m * factorial(l-m) / factorial(l+m) * legendre(l, -m, x)
    else:
        return (2*l-1) / (l+m) * x * legendre(l-1, m, x) - (l-1+m) / (l+m) * legendre(l-2, m, x)

def factorial(n):
    """
    Compute the factorial of n.
    """
    return np.math.factorial(n)