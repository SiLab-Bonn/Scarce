''' Example that shows the potential in a planar sensor in 1D.

    A not fully depleted is supported.The poisson equation is solved
    numerically and compared to the analytical solution.
    The 1D solution is also correct in the 2D case for a sensor with
    100% fill factor.

    The numerical solving for underdepleted sensors is not straight forward,
    due to the a priory unknown depletion depth. This example shows also how
    to deal with that.
'''

import numpy as np
import fipy
import matplotlib.pyplot as plt

from scarce import constant, fields, silicon, solver
from scipy import constants


def calculate_potential(mesh, rho, epsilon, V_read, V_bias, x_dep):
    r''' Calculate the potential with a given space charge distribution.

    If the depletion width is too large the resulting potential will have
    a minimum < bias voltage. This is unphysical.
    '''

    # The field scales with rho / epsilon, thus scale to proper value to
    # counteract numerical instabilities
    epsilon_scaled = 1.
    rho_scale = rho / epsilon

    potential = fipy.CellVariable(mesh=mesh, name='potential', value=0.)

    electrons = fipy.CellVariable(mesh=mesh, name='e-')
    electrons.valence = -1

    electrons.setValue(rho_scale)

    charge = electrons * electrons.valence
    charge.name = "charge"

    # A depletion zone within the bulk requires an internal boundary condition
    # Internal boundary conditions seem to challenge fipy, see:
    # http://www.ctcms.nist.gov/fipy/documentation/USAGE.html#applying-internal-boundary-conditions

    large_value = 1e+15  # Hack for optimizer

    mask = mesh.x > x_dep
    potential.equation = (fipy.DiffusionTerm(coeff=epsilon_scaled) -
                          fipy.ImplicitSourceTerm(mask * large_value) +
                          mask * large_value * V_bias + charge == 0)

    potential.constrain(V_read, mesh.facesLeft)
    potential.constrain(V_bias, mesh.facesRight)

    solver.solve(potential, equation=potential.equation)

    return potential


def get_potential(mesh, rho, epsilon, L, V_read, V_bias, max_iter=10):
    r''' Solves the poisson equation for different depletion depths.

    Iterates until the minimum potential equals the bias potential.
    At this point the depletion width is correctly calculated.
    '''
    x_dep_new = L  # Start with full depletion assumption
    for i in range(max_iter):
        potential = calculate_potential(mesh,
                                        rho=rho,
                                        epsilon=epsilon,
                                        V_read=V_read,
                                        V_bias=V_bias,
                                        x_dep=x_dep_new)

        X = np.array(mesh.getFaceCenters()[0, :])
        x_dep_new = X[np.where(potential == potential.min())][0]

        if (i == 0 and np.allclose(potential.min(), V_bias, rtol=1.e-3)) or \
           (i > 0 and np.allclose(potential.min(), V_bias)):
            return potential

    raise RuntimeError(
        'Depletion region in underdepleted sensor could not be determined')


def create_1D_planar_sensor():
    # Electrical properties
    n_eff = 5e12  # Effective doping concentration
    rho = constants.elementary_charge * n_eff * \
        (1e-4) ** 3  # Charge density in C / um3
    epsilon = constant.epsilon_s * 1e-6  # Permitticity of silicon in F/um

    # External voltages in V
    V_bias = -200.
    V_read = -0

    # Geometry
    dx = 0.01  # Grid spacing / resolution
    L = 200.  # Length of simulation domain / width of sensor in um

    # Create mesh
    nx = L / dx  # Number of space points
    mesh = fipy.Grid1D(dx=np.ones((int(nx), )) * dx, nx=nx)
    X = np.array(mesh.getFaceCenters()[0, :])

    # Get 1D potential with numerical solver
    potential = get_potential(mesh, rho, epsilon, L, V_read, V_bias)

    # Get correct analytical solution
    potential_a = fields.get_potential_planar_analytic_1D(
        X, V_bias=V_bias, V_readout=V_read, n_eff=n_eff, D=L)

    # Get depletion width from empiric formular in um
    x_dep = silicon.get_depletion_depth(
        abs(V_bias), n_eff, temperature=300) * 1e6

    # Plot result
    plt.plot(X, potential_a, '--', linewidth=2, label='Analytic solution')
    plt.plot(np.array(mesh.getFaceCenters()[0, :]), np.array(
        potential.arithmeticFaceValue()), '-', label='Numeric solution')
    plt.plot([x_dep, x_dep], plt.ylim(), '--',
             linewidth=2., label='Depletion zone')
    plt.plot(plt.xlim(), [V_bias, V_bias], '--',
             linewidth=2., label='Bias voltage')
    plt.xlabel('Position [um]')
    plt.ylabel('Potential [V]')
    plt.grid()
    plt.legend(loc=0)
    plt.title('Potential of a %sfully depleted planar silicon sensor\
                with full fill factor' % ('not ' if x_dep < L else ''))
    plt.show()

if __name__ == '__main__':
    create_1D_planar_sensor()
