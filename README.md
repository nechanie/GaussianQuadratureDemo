# GaussianQuadrature

## Overview
This set of Python scripts provides tools for performing numerical integration using various methods including Gauss-Legendre quadrature, the trapezoidal rule, and Simpson's rule. The `GaussQuad.py` file contains classes and functions for calculating basis polynomials and applying Gaussian Quadrature, while `comparisons.py` includes methods for comparing the accuracy of different numerical integration techniques.

---

## `src/GaussQuad.py`

### Description
This file defines classes and methods for computing basis polynomials and performing Gauss Quadrature integration.

### Key Classes and Functions
- `BasisPolynomials`: Class for calculating basis polynomials (Legendre, Chebyshev) and their derivatives.
- `GaussQuadrature`: Class for performing Gauss Quadrature integration.
- `Gaussian_Quad`: Function for executing Gaussian quadrature and returning the integral value, nodes, and weights.

### Usage
1. Import the `GaussQuadrature` class or `Gaussian_Quad` function.
2. Create a function to integrate.
3. Use `Gaussian_Quad` to perform the integration, specifying the number of nodes and the interval.

### CSV File Format
When using `GaussQuadrature.generate_and_save()`, the script generates a CSV file containing two columns:
- **Nodes**: The computed nodes for the quadrature.
- **Weights**: The corresponding weights for each node.

---

## `src/Comparisons.py`

### Description
This file contains functions to compare the performance of various numerical integration methods.

### Key Functions
- `trapezoidal_rule`: Implements the trapezoidal rule for integration.
- `simpsons_rule`: Implements Simpson's rule for integration.
- `compare_methods`: Compares the convergence and error of Gauss-Legendre, trapezoidal rule, and Simpson's rule.

### Usage
1. Import the required functions.
2. Call `compare_methods` with the standard and modified functions, the exact integral value for comparison, the interval of integration, and the maximum power of 2 for the number of intervals or nodes.

### Output
The `compare_methods` function returns the necessary data to create a graph showing the convergence and error of each method.

---

## Requirements
- Python 3.x
- Libraries: `numpy`, `pandas`, `matplotlib`, `scipy`, `math`

### Installation
Ensure that Python 3.x is installed and the required libraries are available. If not, they can be installed via pip:
```bash
pip install numpy pandas matplotlib scipy
