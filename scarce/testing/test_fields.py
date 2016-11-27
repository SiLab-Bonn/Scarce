import unittest
import numpy as np
import os
from fipy import GmshImporter2D
from fipy.tools import dump

from scarce.testing import tools
from scarce import constant
from scarce import fields, geometry


class TestFields(unittest.TestCase):

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

#                 import matplotlib.pyplot as plt
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

        # Create x,y grid
        nx, ny = 1000, 1000
        x = np.linspace(min_x, max_x, nx)
        y = np.linspace(min_y, max_y, ny)
        xx, yy = np.meshgrid(x, y, sparse=True)

        field_description = fields.Description(potential,
                                               min_x=min_x,
                                               max_x=max_x,
                                               min_y=min_y,
                                               max_y=max_y,
                                               nx=nx,
                                               ny=ny,
                                               smoothing=0.2)

        def field_analytic(x, y):
            return fields.get_weighting_field_analytic(x, y, D=thickness, S=width, is_planar=True)

        # Evaluate field on a grid
        f_analytic_x, f_analytic_y = field_analytic(xx, yy)
        f_numeric_x, f_numeric_y = field_description.get_field(xx, yy)

        for i in [-45, -30, -15, -10, 0, 10, 15, 30, 45]:  # Check only at center pixel, edge pixel are not interessting
            self.assertTrue(np.allclose(f_analytic_x.T[nx / 2 + i, :], f_numeric_x.T[nx / 2 + i, :], rtol=0.01, atol=0.01))
            self.assertTrue(np.allclose(f_analytic_y.T[nx / 2 + i, :], f_numeric_y.T[nx / 2 + i, :], rtol=0.01, atol=0.01))
            
    def test_potential_smoothing(self):
        '''  Checks the smoothing of the potential to be independent of the potential values.
        '''

        n_pixel = 11
        width = 50.
        thickness = 50.

        # Create x,y grid
        min_x, max_x = -width * float(n_pixel), width * float(n_pixel)
        min_y, max_y = 0., thickness
        nx, ny = 1000, 1000
        x = np.linspace(min_x, max_x, nx)
        y = np.linspace(min_y, max_y, ny)
        xx, yy = np.meshgrid(x, y, sparse=True)
        
        # Load potential solution to save time
        potential = dump.read(filename = os.path.join(constant.FIXTURE_FOLDER, 'potential.sc'))
        
        def upcale_potential(potential, V_readout, V_bias):
            ''' Scales potential to [V_bias, V_readout] to simulate other bias settings
            '''
            return ((potential - np.nanmin(potential)) / (np.nanmax(potential) - np.nanmin(potential))) * (V_readout - V_bias) + V_readout
        
        def downscale_potential(potential):
            ''' Scales potential to [0, 1] to make the smoothing result comparible
            '''
            return (potential - np.nanmin(potential)) / (np.nanmax(potential) - np.nanmin(potential))

        potential_descr = fields.Description(potential,
                                                   min_x=min_x,
                                                   max_x=max_x,
                                                   min_y=min_y,
                                                   max_y=max_y,
                                                   nx=nx,
                                                   ny=ny)
        
        # Expected result for the std. smoothing value and a potential between 0 and 1
        pot_numeric = downscale_potential(potential_descr.get_potential_smooth(xx, yy))
        
        for V_bias in [-100, -1000]:
            for V_readout in [50, 0, -50]:
                # Create fake data with different bias by upscaling
                potential_scaled = upcale_potential(potential, V_readout, V_bias)
                # Describe upscaled data
                potential_descr_scaled = fields.Description(potential_scaled,
                                                   min_x=min_x,
                                                   max_x=max_x,
                                                   min_y=min_y,
                                                   max_y=max_y,
                                                   nx=nx,
                                                   ny=ny)
                # Downscale smoothed potential for comparison     
                pot_numeric_2 = downscale_potential(potential_descr_scaled.get_potential_smooth(xx, yy))
        
                self.assertTrue(np.allclose(pot_numeric, pot_numeric_2, equal_nan=True))
        

if __name__ == "__main__":
    unittest.main()
