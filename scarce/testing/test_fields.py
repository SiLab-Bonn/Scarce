import unittest
import numpy as np
import os
from fipy import GmshImporter2D

from scarce.testing import tools
from scarce import fields, geometry


class Test(unittest.TestCase):

    @classmethod
    def tearDownClass(cls):
        os.remove('planar_mesh_tmp_2.msh')

    def test_w_potential_analytic(self):
        ''' Check analytic weighting potential calculation
            for planar sensors.
        '''

        self.assertTrue(tools.check_with_fixture(fields.get_weighting_potential_analytic,
                                                 x=np.linspace(-200, 200, 200),
                                                 y=np.linspace(0, 250, 200),
                                                 D=[100, 150, 200, 250],
                                                 S=[5, 50, 100, 150, 250],
                                                 is_planar=True))

    def test_w_field_analytic(self):
        ''' Check analytic weighting field calculation
            for planar sensors.
        '''

        self.assertTrue(tools.check_with_fixture(fields.get_weighting_field_analytic,
                                                 x=np.linspace(-200, 200, 200),
                                                 y=np.linspace(0, 250, 200),
                                                 D=[100, 150, 200, 250],
                                                 S=[5, 50, 100, 150, 250],
                                                 is_planar=True))

    def test_w_pot_field_analytic(self):
        ''' Check weighting potential/field of planar sensor.

        Check if grad(-Phi) = E_w_x, E_w_y
        '''

        for width in [5., 50, 250]:
            for thickness in [50., 100., 200.]:
                x = np.linspace(-width * 2, width * 2, 1000)
                y = np.linspace(0, thickness, 1000)
                xx, yy = np.meshgrid(x, y, sparse=True)

                E_w_x, E_w_y = fields.get_weighting_field_analytic(xx, yy, D=thickness, S=width, is_planar=True)
                Phi_w = fields.get_weighting_potential_analytic(xx, yy, D=thickness, S=width, is_planar=True)

                E_w_y_2, E_w_x_2 = np.gradient(-Phi_w, np.gradient(y), np.gradient(x))

                array_1 = np.array([E_w_x, E_w_y])
                array_2 = np.array([E_w_x_2, E_w_y_2])

                # Assert that less than 1 % of the field poinst have an error > 1%
                self.assertLess(tools.compare_arrays(array_1, array_2, threshold=0.01), 0.01)

    def test_weighting_potential_planar(self):
        '''  Checks the numerical potential estimation by comparing to correct analytical potential.
        '''

        # Influences how correct the field for the center pixel(s) is
        # due to more far away infinite boundary condition
        n_pixel = 11

        for width in [50., 200.]:
            for thickness in [50., 250.]:
                # Analytical solution only existing for pixel width = readout pitch (100% fill factor)
                pitch = width

                print 'width, thickness', width, thickness
                geometry.mesh_planar_sensor(
                    n_pixel=n_pixel,
                    width=width,
                    thickness=thickness,
                    resolution=600. * np.sqrt(width / 50.) * np.sqrt(50. / thickness),
                    filename='planar_mesh_tmp_2.msh')
                mesh = GmshImporter2D('planar_mesh_tmp_2.msh')

                potential = fields.calculate_planar_sensor_w_potential(mesh=mesh,
                                                                       width=width,
                                                                       pitch=pitch,
                                                                       n_pixel=n_pixel,
                                                                       thickness=thickness)

                min_x, max_x = -width * float(n_pixel), width * float(n_pixel)
                min_y, max_y = 0., thickness
                nx, ny = 1000, 1000

                potential_description = fields.Description(potential,
                                                           min_x=min_x,
                                                           max_x=max_x,
                                                           min_y=min_y,
                                                           max_y=max_y,
                                                           nx=nx,
                                                           ny=ny)

                def potential_analytic(x, y):
                    return fields.get_weighting_potential_analytic(x, y, D=thickness, S=width, is_planar=True)

                # Create x,y grid
                x = np.linspace(min_x, max_x, nx)
                y = np.linspace(min_y, max_y, ny)
                xx, yy = np.meshgrid(x, y, sparse=True)

                # Evaluate potential on a grid
                pot_analytic = potential_analytic(xx, yy)
                pot_numeric = potential_description.get_potential(xx, yy)

#                 for i in [0, 10, 15, 30, 45]:
#                     plt.plot(y, pot_analytic.T[nx / 2 + i, :], label='Analytic')
#                 for i in [0, 10, 15, 30, 45]:
#                     plt.plot(y, pot_numeric.T[nx / 2 + i, :], label='Numeric')
#                 plt.legend(loc=0)
#                 plt.show()

                for i in [-45, -30, -15, -10, 0, 10, 15, 30, 45]:  # Check only at center pixel, edge pixel are not interessting
                    sel = pot_analytic.T[nx / 2 + i, :] > 0.01
                    # Check with very tiny and tuned error allowance
                    self.assertTrue(np.allclose(pot_analytic.T[nx / 2 + i, sel], pot_numeric.T[nx / 2 + i, sel], rtol=0.01, atol=0.005))
    @unittest.SkipTest
    def test_weighting_field_planar(self):
        '''  Checks the numerical field estimation by comparing to correct analytical field.
        '''

        width = 50.
        # Analytical solution only existing for pixel width = readout pitch (100 % fill factor)
        pitch = width
        thickness = 200.
        n_pixel = 11

        geometry.mesh_planar_sensor(
            n_pixel=n_pixel,
            width=width,
            thickness=thickness,
            resolution=400,
            filename='planar_mesh_tmp_2.msh')
        mesh = GmshImporter2D('planar_mesh_tmp_2.msh')

        potential = fields.calculate_planar_sensor_w_potential(mesh=mesh,
                                                               width=width,
                                                               pitch=pitch,
                                                               n_pixel=n_pixel,
                                                               thickness=thickness)

        # Define field/potential domain
        min_x, max_x = -width * float(n_pixel), width * float(n_pixel)
        min_y, max_y = 0., thickness

        nx, ny = 1000, 1000
        x = np.linspace(min_x, max_x, nx)
        y = np.linspace(min_y, max_y, ny)

        # Smoothing is needed for gradient forming for E-field
        potential_function = geometry.interpolate_potential(potential, smoothing=0.2)

        # Calculate E-Field from smoothed potential
        field_function_x, field_function_y = geometry.calculate_field(x, y, potential_function)

        def field_analytic(x, y):
            return fields.get_weighting_field_analytic(x, y, D=thickness, S=width, is_planar=True)

        # Create x,y plot grid
        xx, yy = np.meshgrid(x, y, sparse=True)

        # Evaluate potential on a grid
        f_analytic = field_analytic(xx, yy)
        f_numeric_x = field_function_x(xx, yy)
        f_numeric_y = field_function_y(xx, yy)

        for i in [-45, -30, -15, -10, 0, 10, 15, 30, 45]:  # Check only at center pixel, edge pixel are not interessting
            self.assertTrue(np.allclose(f_analytic[0].T[nx / 2 + i, :], f_numeric_x[nx / 2 + i, :], rtol=0.01, atol=0.01))
            self.assertTrue(np.allclose(f_analytic[1].T[nx / 2 + i, :], f_numeric_y[nx / 2 + i, :], rtol=0.01, atol=0.01))

if __name__ == "__main__":
    unittest.main()
