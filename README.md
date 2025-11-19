# python_equation
Algebra Equation Solver (Python)

This project implements a basic algebra equation solver using only the Python standard library. It defines classes for monomials, polynomials, and equations, and supports symbolic representation and simplification. Solving is currently limited to single-variable linear equations.

Usage
from equation import Monomial, Polynomial, Equation

x = Monomial(1, {'x': 1})
p1 = (x + 1) ** 2
p2 = (x + 2) ** 2
eq = Equation(p1, p2)
eq.simplify()


Example output:

-3 - 2x = 0

Project Structure
.
├── equation.py
├── main.py
└── README.md

Reference

Details of implementation are described in the following article:
http://windybird.com/posts/python_equation.html
