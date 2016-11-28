import unittest
import fipy
import numpy as np

from scarce.examples import potential_1D
from scarce import constant, fields
from scipy import constants


class TestSolver(unittest.TestCase):

    def test_linear_poison_solver(self):
        ''' Compare the result of the poison solution with
            analytical result.
        '''

        # Electrical properties
        n_eff = 5e12  # Effective doping concentration
        rho = constants.elementary_charge * n_eff * (1e-4) ** 3  # Charge density in C / um3
        epsilon = constant.epsilon_s * 1e-6  # Permitticity of silicon in F/um

        # External voltages in V
        V_read = -0

        # Geometry
        dx = 0.01  # Grid spacing / resolution
        L = 200.  # Length of simulation domain / width of sensor in um

        # Create mesh
        nx = L / dx  # Number of space points
        mesh = fipy.Grid1D(dx=np.ones((int(nx), )) * dx, nx=nx)
        X = np.array(mesh.getFaceCenters()[0, :])

        for V_bias in range(-10, -200, -20):
            # Get 1D potential with numerical solver
            potential = potential_1D.get_potential(mesh, rho, epsilon, L, V_read, V_bias)

            # Get correct analytical solution
            potential_a = fields.get_potential_planar_analytic_1D(X, V_bias=V_bias, V_readout=V_read, n_eff=n_eff, D=L)

            self.assertTrue(np.allclose(potential, potential_a[:-1], atol=1e-1))

if __name__ == "__main__":
    unittest.main()
