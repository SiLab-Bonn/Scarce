from scarce.constant import epsilon_s

import numpy as np
import time
import fipy

from scarce import constant
from scipy import constants

dx = 0.1
L = 200.  # um

nx = L / dx

# dx = np.linspace(0.01, 0.05, nx)

# The electric fiedl strength scales with rho / permittivity
# Permitivity is a small number and set to 1. here, thus we have to
# change n_eff accordingly. Reason: numerical stability, nothing else

# rho = constants.elementary_charge * 10**12 * (1e-4)**3 # 10^12 / cm3
# 
# epsilon = constant.epsilon_s / 1e-6 # F/um
# 
# print rho / epsilon

n_eff = 0.02
# raise

V_bias = -1.

V_read =  -0

Vbias = None
depletion_depth = L / 1.

mesh = fipy.Grid1D(dx = np.ones((nx,)) * dx, nx = nx)

print mesh.scaledCellVolumes

X = np.array(mesh.getFaceCenters()[0, :])

L = X.max()

# Given the electrostatic potential :math:`\phi`,

potential = fipy.CellVariable(mesh=mesh, name='potential', value=0.)

permittivity = 1.


   
electrons = fipy.CellVariable(mesh=mesh, name='e-')
electrons.valence = -1

x = mesh.cellCenters[0]
electrons.setValue(n_eff / L)
electrons.setValue(0., where=x > depletion_depth)

charge = electrons * electrons.valence
charge.name = "charge"
   
potential.equation = (fipy.DiffusionTerm(coeff = permittivity) + charge == 0)
    
potential.constrain(V_read, mesh.facesLeft)
potential.constrain(V_bias, mesh.facesRight)

# if Vbias is not None:
#     potential.constrain(Vbias, mesh.facesRight)

potential.equation.solve(var=potential)

analytical = fipy.CellVariable(mesh=mesh, name="analytical solution")

def analytic_solution(x):
    
    return n_eff / (2. * permittivity) * (x**2) - n_eff / (permittivity) * L *x
    
    return (x**2)/2 - 2*x

def a_sol(x, V_bias, V_dep, x_dep, D, V_read=0.):
    V = np.zeros_like(x)
    
    const_1 = (V_bias - V_read) / D + x_dep * V_dep / D**2 *(x_dep / D - 2.)
    const_2 = (V_bias - V_read) / D + x_dep**2 * V_dep / D**3
    
    V[x <= x_dep] = V_dep / D**2 * x[x <= x_dep]**2 + const_1 * x[x <= x_dep] + V_read
    V[x > x_dep] = (x[x > x_dep] - D)*const_2 + V_bias
    
    return V

import matplotlib.pyplot as plt

        
  
epsilon = permittivity
x_dep = depletion_depth
        
V_dep = n_eff / (2. * epsilon) * L
print V_dep
        
V = a_sol(X, V_bias=V_bias, V_dep=V_dep, x_dep=x_dep, D=L, V_read=V_read)


plt.plot(X, V, '--', label='analytic')
print mesh.getFaceCenters().shape
plt.plot(np.array(mesh.getFaceCenters()[0, :]), np.array(potential.arithmeticFaceValue()), '-', label='numeric')
plt.legend(loc=0)
plt.show()

# analytical.setValue(analytic_solution(x))

# # Potential is continuous
# analytical.setValue(analytic_solution(x[x.shape[0] / 2 + 1]), where=x > depletion_depth)

# which has been satisifactorily obtained
# 
# print potential.allclose(analytical, rtol = 2e-5, atol = 2e-5)
# 
# viewer = fipy.Viewer(vars=(charge, potential, analytical))#, datamin=0.1, datamax=-1.1)
# viewer.plot()
# raw_input("Press any key to continue...")