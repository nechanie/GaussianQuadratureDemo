# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
import math

# %%
class BasisPolynomials():
    """
    A class for calculating basis polynomials and their derivatives,
    specifically for Legendre polynomials.
    """

    def iterative_legendre(self, n, x):
        """
        Compute Legendre polynomial of degree 'n' at point 'x' iteratively.

        Parameters:
        n (int): Degree of the Legendre polynomial.
        x (float): The point at which the polynomial is evaluated.

        Returns:
        float: The value of the Legendre polynomial of degree 'n' at 'x'.
        """

        P0, P1 = 1, x
        for k in range(2, n + 1):
            Pk = ((2*k - 1)*x*P1 - (k - 1)*P0) / k
            P0, P1 = P1, Pk
        return P0 if n==0 else P1 if n==1 else Pk

    def iterative_legendre_derivative(self, n, x):
        """
        Compute the derivative of the Legendre polynomial of degree 'n' at point 'x' iteratively.

        Parameters:
        n (int): Degree of the Legendre polynomial.
        x (float): The point at which the derivative is evaluated.

        Returns:
        float: The derivative of the Legendre polynomial of degree 'n' at 'x'.
        """

        if n == 0:
            return 0
        else:
            return n * (x * self.iterative_legendre(n, x) - self.iterative_legendre(n - 1, x)) / (x**2 - 1)

    def newton_method(self, n, initial_guess):
        """
        Apply Newton's method to find roots of the Legendre polynomial of degree 'n'.

        Parameters:
        n (int): Degree of the Legendre polynomial.
        initial_guess (float): Initial guess for the root.

        Returns:
        float: An approximate root of the Legendre polynomial of degree 'n'.
        """

        x = initial_guess
        for _ in range(100):
            Pn = self.iterative_legendre(n, x)
            Pn_prime = self.iterative_legendre_derivative(n, x)
            dx = -Pn / Pn_prime
            x += dx
            if abs(dx) < 1e-10:
                break
        return x

    def legendre(self, n, a, b):
        """
        Calculate the nodes and weights for Gauss-Legendre quadrature on interval [a, b].

        Parameters:
        n (int): Number of nodes.
        a (float): Lower bound of the interval.
        b (float): Upper bound of the interval.

        Returns:
        tuple: A tuple containing two lists, the nodes and their corresponding weights.
        """
        nodes = []
        weights = []

        for i in range(1, n + 1):
            # Initial guess for the roots
            initial_guess = math.cos(math.pi * (i - 0.25) / (n + 0.5))

            root = self.newton_method(n, initial_guess)
            
            # Transform the node to the [a, b] interval
            transformed_root = 0.5 * ((b - a) * root + a + b)
            nodes.append(transformed_root)
            
            # Adjust the weight
            weight = 2 / ((1 - root**2) * (self.iterative_legendre_derivative(n, root)**2))
            adjusted_weight = weight * 0.5 * (b - a)
            weights.append(adjusted_weight)

        return nodes, weights



# %%
class GaussQuadrature(BasisPolynomials):
    """
    A class for performing Gauss Quadrature using basis polynomials
    for numerical integration.

    Attributes:
    func (callable): The function to be integrated.
    n (int): Number of nodes.
    interval (tuple): The interval over which to integrate.
    """
    
    def __init__(self, func, n, interval=(-1, 1)):
        """
        Initialize the GaussQuadrature object.

        Parameters:
        func (callable): The function to be integrated.
        n (int): Number of nodes.
        interval (tuple): The interval over which to integrate.
        """
        super().__init__()
        self.func = func
        self.n = n
        self.nodes = None
        self.weights = None
        self.interval = interval

    def generate(self):
        """
        Generate nodes and weights and update class state.
        Returns:
        GaussQuadrature: Returns the instance itself for chaining methods.
        """

        self.nodes, self.weights = self.nodes_and_weights()
        return self
    
    def nodes_and_weights(self):
        """
        Generate nodes and weights based on the chosen method and interval.

        Returns:
        tuple: A tuple containing two lists, the nodes and their corresponding weights.
        """

        return self.legendre(self.n, self.interval[0], self.interval[1])
    
    def gauss_quadrature(self):
        """
        Perform the Gauss quadrature to approximate the integral of the function.

        Returns:
        float: The approximated integral of the function.
        """

        return sum(self.weights * self.func(self.nodes))

def Gaussian_Quad(n, interval, func, filename='nodes_and_weights.csv'):
    """
    Perform Gaussian quadrature for numerical integration.

    Parameters:
    n (int): Number of nodes.
    interval (tuple): The interval over which to integrate.
    func (callable): The function to be integrated.
    filename (str): Filename for saving/loading nodes and weights.

    Returns:
    tuple: Approximated integral value, nodes, and weights.
    """
    
    quadrature = GaussQuadrature(func, n, interval)
    quadrature.generate()
    return quadrature.gauss_quadrature(), quadrature.nodes, quadrature.weights



