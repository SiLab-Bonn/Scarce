from scarce.constant import epsilon_s

import numpy as np
import time
import fipy

from scarce import constant, fields
from scipy import constants

dx = 0.01 # 0.1
L = 1. # 200.  # um

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

n_eff = 1. # 0.02
# raise

V_bias = -0.5

V_read = -0


depletion_depth = L / 1.

mesh = fipy.Grid1D(dx=np.ones((nx,)) * dx, nx=nx)

print mesh.scaledCellVolumes

X = np.array(mesh.getFaceCenters()[0, :])

L = X.max()

# Given the electrostatic potential :math:`\phi`,

potential = fipy.CellVariable(mesh=mesh, name='potential', value=0.)

permittivity = 1.
epsilon = permittivity

electrons = fipy.CellVariable(mesh=mesh, name='e-')
electrons.valence = -1

x = mesh.cellCenters[0]
electrons.setValue(n_eff / L)
electrons.setValue(0., where=x > depletion_depth)

charge = electrons * electrons.valence
charge.name = "charge"

potential.equation = (fipy.DiffusionTerm(coeff=permittivity) + charge == 0)

potential.constrain(V_read, mesh.facesLeft)

x_dep = np.sqrt(2. * epsilon / n_eff * (V_read - V_bias))

if x_dep >= L:
    potential.constrain(V_bias, mesh.facesRight)

potential.equation.solve(var=potential)

analytical = fipy.CellVariable(mesh=mesh, name="analytical solution")


def analytic_solution(x):

    return n_eff / (2. * permittivity) * (x ** 2) - n_eff / (permittivity) * x_dep * x

    return (x ** 2) / 2 - 2 * x


def a_sol(x, V_bias, V_dep, x_dep, D, V_read=0.):
    V = np.zeros_like(x)

    const_1 = (V_bias - V_read) / D + x_dep * V_dep / D ** 2 * (x_dep / D - 2.)
    const_2 = (V_bias - V_read) / D + x_dep ** 2 * V_dep / D ** 3

    V[x <= x_dep] = V_dep / D ** 2 * x[x <= x_dep] ** 2 + const_1 * x[x <= x_dep] + V_read
    V[x > x_dep] = (x[x > x_dep] - D) * const_2 + V_bias

    return V

import matplotlib.pyplot as plt



x_dep = depletion_depth

V_dep = n_eff / (2. * epsilon) * x_dep
print V_dep

V = a_sol(X, V_bias=V_bias, V_dep=V_dep, x_dep=x_dep, D=L, V_read=V_read)

V_2 = fields.get_potential_planar_analytic(X, V_bias=V_bias, V_readout=V_read, epsilon=epsilon, n_eff=n_eff, D=L)
V_3 = analytic_solution(X)
# plt.plot(X, V, '--', label='analytic')
# plt.plot(X, V_3, '--', label='analytic old')
plt.plot(X, V_2, '--', label='analytic new')
print mesh.getFaceCenters().shape
plt.plot(np.array(mesh.getFaceCenters()[0, :]), np.array(potential.arithmeticFaceValue()), '-', label='numeric')
plt.legend(loc=0)
plt.show()

# analytical.setValue(analytic_solution(x))

# Potential is continuous
# analytical.setValue(analytic_solution(x[x.shape[0] / 2 + 1]), where=x > depletion_depth)

# which has been satisifactorily obtained
#
# print potential.allclose(analytical, rtol = 2e-5, atol = 2e-5)
#
# viewer = fipy.Viewer(vars=(charge, potential, analytical))#, datamin=0.1, datamax=-1.1)
# viewer.plot()
# raw_input("Press any key to continue...")