from GaussQuad import Gaussian_Quad
import numpy as np
from scipy.integrate import quad

def trapezoidal_rule(func, a, b, n):
    h = (b - a) / n
    result = 0.5 * (func(a) + func(b))
    for i in range(1, n):
        result += func(a + i * h)
    return result * h

def simpsons_rule(func, a, b, n):
    h = (b - a) / n
    result = func(a) + func(b)
    for i in range(1, n, 2):
        result += 4 * func(a + i * h)
    for i in range(2, n-1, 2):
        result += 2 * func(a + i * h)
    return result * h / 3

def compare_methods(func, interval, max_n):
    gauss_legendre_errors = []
    trapezoidal_errors = []
    simpsons_errors = []
    n_values = [2**k for k in range(max_n + 1)]
    
    exact_integral, _ = quad(func, interval[0], interval[1])
    for n in n_values:
        print(n)
        gauss_legendre_result = Gaussian_Quad(n, interval, func)[0]
        trapezoidal_result = trapezoidal_rule(func, interval[0], interval[1], n)
        simpsons_result = simpsons_rule(func, interval[0], interval[1], n)

        gauss_legendre_errors.append(abs(gauss_legendre_result - exact_integral))
        trapezoidal_errors.append(abs(trapezoidal_result - exact_integral))
        simpsons_errors.append(abs(simpsons_result - exact_integral))
    
    return gauss_legendre_errors, trapezoidal_errors, simpsons_errors

# # Define interval of integration
# a, b = -1, 1

# # Calculate the approximate integral using SciPy's quad function

# # Create visual comparison of various numerical integration methods
# compare_methods(np.exp, (a, b), 7)
