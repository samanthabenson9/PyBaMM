# -*- coding: utf-8 -*-
# from casadi import *
# import numpy as N
# import matplotlib.pyplot as plt


import casadi
import pybamm
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import brentq


"""
Demonstration on how to construct a fixed-step implicit Runge-Kutta integrator
@author: Joel Andersson, K.U. Leuven 2013
"""

# End time
tf = 1.0

# Dimensions
n_x = 1
n_p = 0

# Declare variables
x = casadi.SX.sym("x", n_x)  # state
p = casadi.SX.sym("u", n_p)  # control

# ODE right hand side function
# ode = vertcat((1 - x[1]*x[1])*x[0] - x[1] + p, \
#   x[0], \
#   x[0]*x[0] + x[1]*x[1] + p*p)

ode = 0.1 * x[0]
dae = {"x": x, "p": p, "ode": ode}
f = casadi.Function("f", [x, p], [ode])

# Number of finite elements
n = 100

# Size of the finite elements
h = tf / n

# Degree of interpolating polynomial
d = 4

# Choose collocation points
tau_root = [0] + casadi.collocation_points(d, "legendre")

# Coefficients of the collocation equation
C = np.zeros((d + 1, d + 1))

# Coefficients of the continuity equation
D = np.zeros(d + 1)

# Dimensionless time inside one control interval
tau = casadi.SX.sym("tau")

# For all collocation points
for j in range(d + 1):
    # Construct Lagrange polynomials to get the polynomial basis at the collocation point
    L = 1
    for r in range(d + 1):
        if r != j:
            L *= (tau - tau_root[r]) / (tau_root[j] - tau_root[r])

    # Evaluate the polynomial at the final time to get the coefficients of the continuity equation
    lfcn = casadi.Function("lfcn", [tau], [L])
    D[j] = lfcn(1.0)

    # Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
    tfcn = casadi.Function("tfcn", [tau], [casadi.tangent(L, tau)])
    for r in range(d + 1):
        C[j, r] = tfcn(tau_root[r])

# Total number of variables for one finite element
X0 = casadi.MX.sym("X0", n_x)
P = casadi.MX.sym("P", n_p)
V = casadi.MX.sym("V", d * n_x)

# Get the state at each collocation point
X = [X0] + casadi.vertsplit(V, [r * n_x for r in range(d + 1)])

# Get the collocation quations (that define V)
V_eq = []
for j in range(1, d + 1):
    # Expression for the state derivative at the collocation point
    xp_j = 0
    for r in range(d + 1):
        xp_j += C[r, j] * X[r]

    # Append collocation equations
    f_j = f(X[j], P)
    V_eq.append(h * f_j - xp_j)

# Concatenate constraints
V_eq = casadi.vertcat(*V_eq)

# Root-finding function, implicitly defines V as a function of X0 and P
vfcn = casadi.Function("vfcn", [V, X0, P], [V_eq])

# Convert to SX to decrease overhead
vfcn_sx = vfcn.expand()

# Create a implicit function instance to solve the system of equations
ifcn = casadi.rootfinder("ifcn", "fast_newton", vfcn_sx)
V = ifcn(casadi.MX(), X0, P)
X = [X0 if r == 0 else V[(r - 1) * n_x : r * n_x] for r in range(d + 1)]

# Get an expression for the state at the end of the finie element
XF = 0
for r in range(d + 1):
    XF += D[r] * X[r]

# Get the discrete time dynamics
F = casadi.Function("F", [X0, P], [XF])

# Do this iteratively for all finite elements
Xs = casadi.MX(n_x, n+1)
Xs[:, 0] = X0
X = X0
for i in range(0,n):
    X = F(X, P)
    XZ = F(Xs[:, i], P)
    Xs[:, i + 1] = XZ[0]
# X = F(X, P)

# Fixed-step integrator
irk_integrator = casadi.Function(
    "irk_integrator", {"x0": X0, "p": P, "xf": Xs}, casadi.integrator_in(), casadi.integrator_out()
)
sol_irk =irk_integrator(x0=1, p=[])
print(sol_irk['xf'])