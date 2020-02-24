#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 09:07:07 2019

@author: duncanmartinson
"""
## You must download Docker to run this program
## Run these lines of code once you have downloaded docker:

# cd (directory)
# docker run -ti -v $(pwd):/home/fenics/shared:z quay.io/fenicsproject/stable
# cd $HOME/shared
# python3 Pillay_2D_Model.py
###############################################################################
# Solves the 2D Pillay PDE System described in Pillay et al. (2017)
# on a Unit Square domain with No-Flux BCs and either (1) several
# Delta Dirac Point sources as the initial condition, or (2) data from a CSV 
# file containing the results of a 2D CA model.
# Parameters are taken from Pillay et al. (2017), as well as Byrne and Chaplain
# (1995)
###############################################################################
from dolfin import *
from mshr import *
from math import *
import numpy as np
from scipy.interpolate import NearestNDInterpolator
from decimal import Decimal

# Read csv data for tip and stalk cells
TCdata = np.loadtxt('10feb2020_PillayCAModel_Fixed_P_p0_P_m1_k100_LinearTAFField_NoBranchingOrAnastomosis_1000Realizations_Average2DTipCellNetwork.csv', delimiter = ',')
ECdata = np.loadtxt('10feb2020_PillayCAModel_Fixed_P_p0_P_m1_k100_LinearTAFField_NoBranchingOrAnastomosis_1000Realizations_Average2DStalkCellNetwork.csv', delimiter = ',')

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
# using the .csv file obtained by averaging over realizations of the CA Model:
class InitialConditionfromCSVFile(UserExpression):
    def __init__(self, f, **kwargs):
        self._f = f
        UserExpression.__init__(self, **kwargs)
    def eval(self, values, x):
        values[:] = self._f(*x)

# If c(x,y,t) is not included in the variational statement of the problem, 
# these lines of code allow you to formulate it so as to save 
# computational effort (this only works if c(x,y,t) is at steady state).
# Alternatively, it can be used to set up the initial condition for c
class TAF_Profile(UserExpression):
    def eval(self, values, x):
        #values[0] = x[0] # c(x,y,t) = x
        values[:] = 0.5*(x[0] + x[1]) # c(x,y,t) = 0.5*(x+y)

# Number of intervals in x and y direction:
n_x = 200
n_y = 200
# Define domain and mesh:
# Unit Square Mesh
mesh = UnitSquareMesh(n_x,n_y)
(x,y) = SpatialCoordinate(mesh) # Define labels for Cartesian coordinates

# Define Finite Elements for Functions:
# Quadratic Lagrange polynomials on triangle element
V = FiniteElement("CG", triangle, 2) # Define element
# Create a mixed finite element to solve for n(x,y,t), e(x,y,t), c(x,y,t)
# simultaneously:
A = FunctionSpace(mesh, V) # This is a single scalar function space for c
Q = FunctionSpace(mesh, MixedElement([V,V]))

# Define a function q which stores n(x,y,t), e(x,y,t), c(x,y,t)
q = Function(Q)

# Create a function incorporating the values of q_{t-1}
q_prevs = Function(Q)
# Define an arbitrary test function q on the mixed Function Space:
p = TestFunction(Q)
# Split this test function into 3 separate ones for n(x,y,t), e(x,y,t), c(x,y,t)
(v,r) = split(p) # Define arbitrary test functions

# Set up Initial Condition on the domain:
# Using the CSV Initial Condition:
print('Interpolating initial condition onto FEM mesh')
print('Tip cells...')
n_init_cond = interpolate(InitialConditionfromCSVFile(TCinterpolant, element = A.ufl_element()), A)
print('Done!')
print('Stalk Cells...')
e_init_cond = interpolate(InitialConditionfromCSVFile(ECinterpolant, element = A.ufl_element()), A)
print('Done!')


# Set up Initial Condition on the domain:
# Using a User-specified initial condition:
assign(q, [n_init_cond, e_init_cond])
assign(q_prevs, q)

print('TAF Concentration...')
## Set up c(x,y,t). Since this is assumed to be at steady state (c(x,y,t) = x),
## this is not included in the mixed Finite Element
c = Function(A)
c.assign(interpolate(TAF_Profile(element = A.ufl_element()), A))
print('Done!')

# Split the mixed Function into 3 separate functions, where u = n, w = e, c = c:
(n,e) = split(q)

# Define a function which evaluates the inner product of the gradients of two scalar fields:
def rhs(a,b):
    return -inner(grad(a), grad(b))*dx

# Enter all parameters as positive numbers 
# (negative values already have their signs considered in the variational form below)
D = Constant(0.001) # Diffusion coefficient for tip cells
chi = Constant(0.4) # Chemotactic sensitivity for tip cells
Lambda = Constant(0.0) # Branching Rate of tip cells
mu = Constant(160.0) # Rate of Anastomosis, TC Movement (Diffusion and Chemotaxis combined)
a_n = Constant(0.0) # Tip-to-Tip Anastomosis is active
a_e = Constant(0.0) # Tip-to-Sprout Anastomosis is active but set to be a low value.
D_b = Constant(0.0) # y-direction branching rate


# Set up time values:
T = Decimal("2.0") # Final Time
t = Decimal("0.2") # Initial Time
h = Decimal("0.001") # Time Step
dt = Constant(float(h)) # Time Step Represented as a floating point number
nTimeStep = 0 # Number of time steps taken so far

# Modified Snail-Trail Model as Presented in Pillay et al. (2017),
# Crank-Nicolson Variational Formulation:

F_n = (n*v*dx-q_prevs[0]*v*dx
        -dt*D*rhs(0.5*(n+q_prevs[0]),v) 
        - dt*chi*(0.5*(n+q_prevs[0])*inner(grad(c), grad(v)))*dx
        - dt*Lambda*0.5*(n+q_prevs[0])*0.5*c*v*dx
        + dt*mu*a_n*0.5*(n+q_prevs[0])*0.5*(n+q_prevs[0])*v*dx
        + dt*mu*a_e*0.5*(n+q_prevs[0])*0.5*(e+q_prevs[1])*v*dx
        - dt*D_b*grad(0.5*(n+q_prevs[0])*grad(c)[1])[1]*v*dx
        - dt*D_b*grad(c*grad(0.5*(n+q_prevs[0]))[1])[1]*v*dx
        + dt*(a_n*0.5*(n+q_prevs[0]) + a_e*0.5*(e+q_prevs[1]))*div(D*grad(0.5*(n+q_prevs[0]))-chi*0.5*(n+q_prevs[0])*grad(c))*v*dx
        )
F_e = (e*r*dx -q_prevs[1]*r*dx 
        -dt*mu*0.5*(n+q_prevs[0])*r*dx - dt*mu*a_n*0.25*(n+q_prevs[0])*(n+q_prevs[0])*r*dx + dt*a_n*0.5*(n+q_prevs[0])*div(D*grad(0.5*(n+q_prevs[0]))-chi*0.5*(n+q_prevs[0])*grad(c))*r*dx
      )

F_cn = F_n + F_e
# Create output .pvd Files to store the solution of n(x,y,t), e(x,y,t), 
# c(x,y,t) at various time values:
output1 = File("10feb2020_PillayPDEModel_Fixed_P_p0_P_m1_k100_TAFFieldLinearinXandY_NoBranchingOrAnastomosis_1000Realizations_Tip_Cells_CSVInitCond.pvd") # n(x,y,t)
output2 = File("10feb2020_PillayPDEModel_Fixed_P_p0_P_m1_k100_TAFFieldLinearinXandY_NoBranchingOrAnastomosis_1000Realizations_Stalk_Cells_CSVInitCond.pvd") # e(x,y,t)
output3 = File("10feb2020_PillayPDEModel_Fixed_P_p0_P_m1_k100_TAFFieldLinearinXandY_NoBranchingOrAnastomosis_1000Realizations_TAF_Concentration_CSVInitCond.pvd") # c(x,y,t)

# Output initial condition to files
(n0,e0) = q_prevs.split()
output1<<(n0, float(t))
output2<<(e0, float(t))
output3<<(c, float(t))


while True:
    # Update time step
    t += h
    print("Solving for Time: ", float(t))
###########################   Crank-Nicolson Method      ######################
    # Run Crank-Nicolson step
    solve(F_cn==0, q)
    
    # Reset q_{n-1}:
    q_prevs.assign(q)
    
    # Add up time step
    nTimeStep += 1
    
    # Output to files only in increments of 0.01 seconds
    if nTimeStep%10==0:
        (n_k, e_k) = q_prevs.split()
        output1<<(n_k, float(t))
        output2<<(e_k, float(t))
    
    if t >= T: break