#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 13:05:17 2019

@author: duncanmartinson
"""
## You must download Docker to run this program
## Run these lines of code once you have downloaded docker and fenics:

# cd (directory)
# docker run -ti -v $(pwd):/home/fenics/shared:z quay.io/fenicsproject/stable
# cd $HOME/shared
# python3 SnailTrail_2D_PDE.py
###############################################################################
# Solves the Snail-Trail PDE System described in Martinson et al. (2020)
# in 2D on a Unit Square domain with No-Flux BCs and data from a CSV
# file containing the results of a 2D CA model as the initial condition.
###############################################################################
from dolfin import *
from mshr import *
from math import exp
import numpy as np
from scipy.interpolate import NearestNDInterpolator
from decimal import Decimal

# Read csv data for tip and stalk cells
TCdata = np.loadtxt('', delimiter = ',')
ECdata = np.loadtxt('', delimiter = ',')

# Store the data you have read in arrays, and use this to create an interpolant
# function on a 2D grid
X, Y, TCValues = TCdata[:,0], TCdata[:,1], TCdata[:,2]
ECValues = ECdata[:,2]
print('Finding interpolant function for TC data...')
TCinterpolant = NearestNDInterpolator((X,Y), TCValues)
print('Done!')
print('Finding interpolant function for EC data...')
ECinterpolant = NearestNDInterpolator((X,Y), ECValues)
print('Done!')
# This class allows the user to set up an initial condition for the PDE model
# using the .csv file:
class InitialConditionfromCSVFile(UserExpression):
    def __init__(self, f, **kwargs):
        self._f = f #Define the interpolator as part of your object
        UserExpression.__init__(self, **kwargs)
    def eval(self, values, x):
        values[:] = self._f(*x)+(10**(-16)) # Interpolate any value of (x,y) on the domain using the interpolator, and perturb values slightly by floating point error to ensure well-conditioned nonlinear solves

# If c(x,y,t) is not included in the variational statement of the problem,
# these lines of code allow you to formulate it so as to save
# computational effort (this only works if c(x,y,t) is at steady state).
# Alternatively, it can be used to set up the initial condition for c
class TAF_Profile(UserExpression):
    def eval(self, values, x):
        values[:] = x[0] # c(x,y,t) = x
        #values[:] = exp(-(x[0]-(1+10**(-16)))**2-(x[1]-1.0)**2/2.5)
        #values[:] = 0.5*(x[1]) + 1.0/3.0*(x[0])
        #values[:] = 0.5*(x[0]+x[1])
        #values[:] = x[0]*x[1]
        #values[:] = 1-(x[0]-0.5)**2-(x[1]-0.5)**2
        #values[:] = 1.0

class TipCell_Initial_Condition(UserExpression):
    def eval(self, values, x):
        values[:] = exp(-x[0]**2/5.0/pow(10,-3))*sin(6*pi*x[1])**2

class StalkCell_Initial_Condition(UserExpression):
    def eval(self, values, x):
        values[:] = 0.0

# Number of finite elements in x and y direction:
n_x = 200
n_y = 200

# Define domain and mesh:
mesh = UnitSquareMesh(n_x,n_y)
# Define labels for Cartesian coordinates
(x,y) = SpatialCoordinate(mesh)

# Define Finite Elements for Functions:
# Quadratic Lagrange polynomials on triangle element
V = FiniteElement("CG", triangle, 2)

# Create a mixed finite element to solve for n(x,y,t) and e(x,y,t)
# simultaneously:
A = FunctionSpace(mesh, V) # This is a single scalar function space for either n, e, or c
Q = FunctionSpace(mesh, MixedElement([V,V]))

# Define a function q which stores [n(x,y,t), e(x,y,t)]
q = Function(Q)

# Create a function incorporating the values of q_{t-1}
q_prevs = Function(Q)

# Define an arbitrary test function on the mixed Function Space:
p = TestFunction(Q)
# Split this test function into 2 separate ones for n(x,y,t), e(x,y,t):
(v,r) = split(p)

# Set up Initial Condition on the domain:
print('Interpolating initial condition onto FEM mesh')
print('Tip cells...')
# Uncomment the following line of code if you wish to specify the initial condition using ABM data
n_init_cond = interpolate(InitialConditionfromCSVFile(TCinterpolant, element = A.ufl_element()), A)
# Uncomment the following line of code if you wish to specify a custom initial condition
#n_init_cond = interpolate(TipCell_Initial_Condition(element=A.ufl_element()),A)
print('Done!')
print('Stalk Cells...')
# Uncomment this line of code if you wish to specify the initial condition using ABM data
e_init_cond = interpolate(InitialConditionfromCSVFile(ECinterpolant, element = A.ufl_element()), A)
# Uncomment the following line of code if you wish to specify a custom initial condition
#e_init_cond = interpolate(StalkCell_Initial_Condition(element=A.ufl_element()),A)
print('Done!')

# Assign to each component in q and q_prevs the functions you have just interpolated
assign(q, [n_init_cond, e_init_cond])
assign(q_prevs, [n_init_cond, e_init_cond])

## Set up c(x,y). Since this is assumed to be at steady state, it is
## not included in the mixed Finite Element
print('TAF Concentration...')
c = Function(A)
c.assign(interpolate(TAF_Profile(element = A.ufl_element()), A))
print('Done!')

# Split the mixed Function q into 2 separate functions: n is for tip cells, and e is for stalk cells
(n,e) = split(q)

# Define a function which evaluates the inner product of the gradients of two scalar fields:
def rhs(a,b):
    return -inner(grad(a), grad(b))*dx

# Enter all parameters as positive numbers
# (negative values already have their signs considered in the variational form below)
D = Constant(0.001) # Random motility coefficient for tip cells
chi = Constant(0.4) # Chemotactic sensitivity for tip cells
Lambda = Constant(0.16) # Branching Rate of tip cells
beta_n = Constant(160.0) # Rate of Tip-to-Tip Anastomosis
beta_e = Constant(5.0) # Rate of Tip-to-Sprout Anastomosis
Dx = Constant(1.0/200) # Spatial Step Size of Agent-Based Model/Length of a Cell
Kappa = 4.0*D/chi/sqrt(inner(grad(c),grad(c)))/Dx
#Kappa = Constant(2.0)

# Uncomment the following line of code for cases in which Kappa is a constant,
# to check that you have the right value:
#print("Kappa = ", float(Kappa))

# Set up time values:
T = Decimal("2.0") # Final Time
t = Decimal("0.0") # Initial Time
h = Decimal("0.001") # Time Step
dt = Constant(float(h)) # Time Step Represented as a floating point number
nTimeStep = 0 # Number of time steps taken so far

dxm = dx()

# Modified Snail-Trail Model,
# Crank-Nicolson Variational Formulation:
F_cn = (n*v*dxm -q_prevs[0]*v*dxm
        -dt*D*rhs(0.5*(n+q_prevs[0]),v)
        - dt*chi*(0.5*(n+q_prevs[0])*inner(grad(c), grad(v)))*dxm
        - dt*Lambda*0.5*(n+q_prevs[0])*c*v*dxm
        + dt*beta_n*0.5*(n+q_prevs[0])*0.5*(n+q_prevs[0])*v*dxm
        + dt*beta_e*0.5*(n+q_prevs[0])*0.5*(e+q_prevs[1])*v*dxm
        +(e-q_prevs[1])*r*dxm
        - dt*Kappa/Dx*sqrt(inner(D*grad(0.5*(n+q_prevs[0]))- chi*0.5*(n+q_prevs[0])*grad(c), D*grad(0.5*(n+q_prevs[0]))- chi*0.5*(n+q_prevs[0])*grad(c)))*r*dxm
        )

# Create output .pvd Files to store the solution of n(x,y,t), e(x,y,t),
# c(x,y,t) at various time points:

output1 = File("Tip_Cells.pvd") # n(x,y,t)
output2 = File("Stalk_Cells.pvd") # e(x,y,t)
output3 = File("TAF_Concentration.pvd")

# Output initial condition to files
(n0,e0) = q_prevs.split()
output1<<(n0, float(t))
output2<<(e0, float(t))
output3<<(c, float(t))


# Solve the PDE system using Finite Elements:
while True:
    # Update time step
    t += h
    print("Solving for Time: ", float(t))
###########################   Crank-Nicolson Method      ######################
    # Solve the nonlinear PDE system using Crank-Nicolson time-stepping and GMRES for linear solver, with Jacobi preconditioner:
    solve(F_cn==0, q, solver_parameters={'newton_solver': {'linear_solver':'gmres', 'preconditioner': 'jacobi'}})

    # Reset q_{t-1}
    q_prevs.assign(q)

    # Increase time step counter
    nTimeStep += 1

    # Output to files only in increments of 0.01 seconds
    if nTimeStep%10==0:
        (n_k, e_k) = q_prevs.split()
        output1<<(n_k, float(t))
        output2<<(e_k, float(t))
    if t >= T: break
