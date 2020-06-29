## You must download Docker to run this program
## Run these lines of code once you have downloaded docker and fenics:

# cd (directory)
# docker run -ti -v $(pwd):/home/fenics/shared:z quay.io/fenicsproject/stable
# cd $HOME/shared
# python3 SnailTrail_L2Norm_Optimize_Beta_e.py
###############################################################################
# Solves for the optimal value of Beta_e used in the 2D Snail-Trail PDE System
# described in Martinson et al. (2020) for a Unit Square domain with No-Flux
# BCs.
# The value of Beta_e minimizes 2D ABM data specified as CSV files. ABM Data at
# t = 0.2 is used as an initial condition for the PDE system, which is solved
# for t = 0.2 to t = 2.
###############################################################################
from dolfin import *
from mshr import *
from math import exp
import numpy as np
from scipy.interpolate import NearestNDInterpolator
from scipy.optimize import least_squares
from decimal import Decimal

# Read csv data for tip and stalk cells at t=0.2s
TCdata = np.loadtxt('Tip_Cell_ABM_Data_t2.csv', delimiter = ',')
ECdata = np.loadtxt('Stalk_Cell_ABM_Data_t2.csv', delimiter = ',')

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

def y_fun(*args, **kwargs):
    TCdata = np.loadtxt('Tip_Cell_ABM_Data_t2.csv', delimiter = ',')
    ECdata = np.loadtxt('Stalk_Cell_ABM_Data_t2.csv', delimiter = ',')
    TCdata = TCdata[:,2]
    ECdata = ECdata[:,2]
    TCdata = np.reshape(np.transpose(np.reshape(TCdata, (201,201))), 201**2)
    TCdata = np.reshape(np.transpose(np.reshape(ECdata, (201,201))), 201**2)
    for i in [4,6,8,10,12,14,16,18,20]:
        TC_name = 'Tip_Cell_ABM_Data_t'+str(i)+'.csv'
        EC_name = 'Stalk_Cell_ABM_Data_t'+str(i)+'.csv'
        next_TCdata = np.loadtxt(TC_name, delimiter = ',')
        next_TCdata = next_TCdata[:,2]
        next_ECdata = np.loadtxt(EC_name, delimiter = ',')
        next_ECdata = next_ECdata[:,2]
        next_TCdata = np.reshape(np.transpose(np.reshape(next_TCdata, (201,201))), 201**2)
        next_ECdata = np.reshape(np.transpose(np.reshape(next_ECdata, (201,201))), 201**2)
        TCdata = np.concatenate((TCdata,next_TCdata))
        ECdata = np.concatenate((ECdata,next_ECdata))
    return (TCdata, ECdata)

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
        #values[:] = exp(-(x[0]-(1.0))**2-(x[1]-1.0)**2/2.5)
        #values[:] = 0.5*(x[1]) + 1.0/3.0*(x[0])
        #values[:] = 0.5*(x[0]+x[1])
        #values[:] = x[0]*x[1]
        #values[:] = 1-(x[0]-0.5)**2-(x[1]-0.5)**2
        #values[:] = 1.0

# Number of finite elements in x and y direction:
n_x = 200
n_y = 200

# Define the domain
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

# Set up Initial Condition on the domain:
# Using the CSV Initial Condition:
print('Interpolating initial condition onto FEM mesh')
print('Tip cells...')
n_init_cond = interpolate(InitialConditionfromCSVFile(TCinterpolant, element = A.ufl_element()), A)
print('Done!')
print('Stalk Cells...')
e_init_cond = interpolate(InitialConditionfromCSVFile(ECinterpolant, element = A.ufl_element()), A)
print('Done!')

## Set up c(x,y,t). Since this is assumed to be at steady state, it is
## not included in the mixed Finite Element
print('TAF Concentration...')
c = Function(A)
c.assign(interpolate(TAF_Profile(element = A.ufl_element()), A))
print('Done!')

# Record ABM Data at time points t = 0.2, 0.4, ..., 2 using the function
# defined above:
ABM_TCdata,ABM_ECdata = y_fun()

# Define a function which evaluates the inner product of the gradients of two scalar fields:
def rhs(a,b):
    return -inner(grad(a), grad(b))*dx

# Define a function that calculates the residual between the ABM and PDE results
def res_fun_L2norm(Beta_e, ABM_TCdata, ABM_ECdata, mesh, x, y, V, A, Q, n_init_cond, e_init_cond, c):
    # Define a function q which stores n(x,y,t), e(x,y,t)
    q = Function(Q)

    # Create a function incorporating the values of q_{t-1}
    q_prevs = Function(Q)

    # Define an arbitrary test function on the mixed Function Space:
    p = TestFunction(Q)
    # Split this test function into 2 separate ones for n(x,y,t), e(x,y,t):
    (v,r) = split(p) # Define arbitrary test functions


    # Assign to each component in q and q_prevs the functions you have just interpolated
    assign(q, [n_init_cond, e_init_cond])
    assign(q_prevs, [n_init_cond, e_init_cond])

    # Split the mixed Function q into 2 separate functions: n is for tip cells, and e is for stalk cells
    (n,e) = split(q)

    # Enter all parameters as positive numbers
    # (negative values already have their signs considered in the variational form below)
    D = Constant(0.001) # Random motility coefficient for tip cells
    chi = Constant(0.4) # Chemotactic sensitivity for tip cells
    Lambda = Constant(0.16) # Branching Rate of tip cells
    beta_n = Constant(160.0) # Rate of Tip-to-Tip Anastomosis
    beta_e = Constant(Beta_e[0]) # Rate of Tip-to-Sprout Anastomosis
    Dx = Constant(1.0/200) # Spatial Step Size of Agent-Based Model
    Kappa = Constant(2.0)

    # Print the current value of beta_e to the screen:
    print('Beta_e = ', float(beta_e))

    # Uncomment the following line of code for cases in which Kappa is a constant,
    # to check that you have the right value:
    # print("Kappa = ", float(Kappa))

    # Set up time values:
    T = Decimal("2.0") # Final Time
    t = Decimal("0.2") # Initial Time
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
    output1 = File("Tip_Cells_Optimized_Beta_e.pvd") # n(x,y,t)
    output2 = File("Stalk_Cells_Optimized_Beta_e.pvd") # e(x,y,t)
    output3 = File("TAF_Concentration_Optimized_Beta_e.pvd") # c(x,y,t)

    # Output initial condition to files
    (n0,e0) = q_prevs.split()
    output1<<(n0, float(t))
    output2<<(e0, float(t))
    output3<<(c, float(t))

    TC_CSVname = 'Tip_Cells_Optimized_Beta_e_CSV_Data_t'
    EC_CSVname = 'Stalk_Cells_Optimized_Beta_e_CSV_Data_t'
    X2 = mesh.coordinates()
    save_TC_name = TC_CSVname + str(2) + '.csv'
    save_EC_name = EC_CSVname + str(2) + '.csv'
    np.savetxt(save_TC_name, (X2[:,0], X2[:,1], n0.compute_vertex_values(mesh)), delimiter=',')
    np.savetxt(save_EC_name, (X2[:,0], X2[:,1], e0.compute_vertex_values(mesh)), delimiter=',')

    PDE_TCdata = n0.compute_vertex_values(mesh)
    PDE_ECdata = e0.compute_vertex_values(mesh)
    # Solve the PDE system using Finite Elements:
    while True:
        # Update time step
        t += h
        print("Solving for Time: ", float(t))
###########################   Crank-Nicolson Method      ######################
        # Solve the nonlinear PDE system using Crank-Nicolson time-stepping and GMRES for linear solver, with Jacobi preconditioner:
        solve(F_cn==0, q, solver_parameters={'newton_solver': {'linear_solver':'gmres', 'preconditioner': 'jacobi'}})

        # Reset q_{n-1}
        q_prevs.assign(q)

        # Increase time step counter
        nTimeStep += 1

        # Output to PVD files in increments of 0.01 seconds
        if nTimeStep%10==0:
            (n_k, e_k) = q_prevs.split()
            output1<<(n_k, float(t))
            output2<<(e_k, float(t))
        # Output to CSV files in increments of 0.2 seconds
        if nTimeStep%200==0:
            (n_k, e_k) = q_prevs.split()
            i = int(round(nTimeStep/200)*2+2)
            # Save solution
            save_TC_name = TC_CSVname + str(i) + '.csv'
            save_EC_name = EC_CSVname + str(i) + '.csv'
            np.savetxt(save_TC_name, (X2[:,0], X2[:,1], n_k.compute_vertex_values(mesh)), delimiter=',')
            np.savetxt(save_EC_name, (X2[:,0], X2[:,1], e_k.compute_vertex_values(mesh)), delimiter=',')
            PDE_TCdata = np.concatenate((PDE_TCdata, n_k.compute_vertex_values(mesh)))
            PDE_ECdata = np.concatenate((PDE_ECdata, e_k.compute_vertex_values(mesh)))
        if t >= T: break

    # Once the loop is broken, calculate the residuals:
    TC_res = ABM_TCdata-PDE_TCdata
    EC_res = ABM_ECdata-PDE_ECdata
    return np.concatenate((TC_res,EC_res))

# Initial guess for beta_e:
beta_e0 = 5.0
# Solve for beta_e using a nonlinear least squares method:
optimal_beta_e = least_squares(res_fun_L2norm, beta_e0, bounds = (0.0,160.0), args = (ABM_TCdata, ABM_ECdata, mesh, x, y, V, A, Q, n_init_cond, e_init_cond, c))

# Print to terminal the value of beta_e and the termination conditions
# of the algorithm
print('Optimal value of Beta_e is ')
print(optimal_beta_e.x)
print(optimal_beta_e.message)

# Save the value of beta_e, the residuals for that value, and the jacobian
# so that you can compute 95% confidence intervals:
np.savetxt('Optimal_Beta_e_Value.csv', optimal_beta_e.x, delimiter=',')
np.savetxt('Optimal_Beta_e_Residuals.csv', optimal_beta_e.fun, delimiter=',')
np.savetxt('Optimal_Beta_e_Jacobian.csv', optimal_beta_e.jac, delimiter=',')
