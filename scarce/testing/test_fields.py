import unittest
import numpy as np
from fipy import GmshImporter2D

from scarce.testing import tools
from scarce import fields, geometry, plot

import matplotlib.pyplot as plt
from matplotlib import colors, cm

from scipy.interpolate import RectBivariateSpline, SmoothBivariateSpline, interp2d

class Test(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

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
#                 y = np.square(np.linspace(0, np.sqrt(thickness), 1000))
                xx, yy = np.meshgrid(x, y, sparse=True)
        
                E_w_x, E_w_y = fields.get_weighting_field_analytic(xx, yy, D=thickness, S=width, is_planar=True)
                Phi_w = fields.get_weighting_potential_analytic(xx, yy, D=thickness, S=width, is_planar=True)
        
                E_w_x_2, E_w_y_2  = np.gradient(-Phi_w, np.gradient(x), np.gradient(y), axis=(1, 0))
                
                array_1 = np.array([E_w_x, E_w_y])
                array_2 = np.array([E_w_x_2, E_w_y_2])

                # Assert that less than 1 % of the field poinst have an error > 1%
                self.assertLess(tools.compare_arrays(array_1, array_2, threshold=0.01), 0.01) 
    
    @unittest.SkipTest            
    def test_weighting_field_planar(self):
        '''  Checks the numerical field estimation by comparing to correct analytical field.
        '''
        
        width = 50.
        # Analytical solution only existing for pixel width = readout pitch (100 % fill factor)
        pitch = width
        thickness = 200.
        n_pixel = 9

#         mesh = geometry.mesh_planar_sensor(
#                              n_pixel=n_pixel,
#                              width=width, 
#                              thickness=thickness, 
#                              resolution=200)
         
#         plot.plot_mesh(mesh)
#         raise
        mesh = GmshImporter2D('sensor.msh')
        
        potential = fields.calculate_planar_sensor_w_potential(mesh=mesh, 
                                                 width=width, 
                                                 pitch=pitch,
                                                 n_pixel=n_pixel, 
                                                 thickness=thickness)

        potential_function = geometry.interpolate_potential(potential)
        
        def potential_analytic(x, y):
            return fields.get_weighting_potential_analytic(x, y, D=thickness, S=width, is_planar=True)
        
        min_x, max_x=-width * float(n_pixel) / 2., width * float(n_pixel) / 2.
        min_y, max_y = 0., thickness
        
        nx, ny = 1000, 1000
        x = np.linspace(min_x, max_x, nx)
        y = np.linspace(min_y, max_y, ny)
 
        # Create x,y plot grid
        xx, yy = np.meshgrid(x, y)#, sparse=True)
        
#         print potential_function(xx, yy)
#         raise
        
#         potential_function_smooth = RectBivariateSpline(y, x, potential_function(xx, yy).T)

        zz = np.array(potential_function(xx, yy)).ravel()
#         print z.shape
#         print type(z)
#         z=np.zeros((x.shape[0], y.shape[0]))
#         print x.shape, y.shape, z.shape

        sel = zz == np.isfinite

        print zz[sel]
        raise
        potential_function_smooth = SmoothBivariateSpline(np.array(xx.ravel())[sel], np.array(yy.ravel())[sel], np.array(zz.ravel())[sel], kx=3, ky=3)
        
        print potential_function_smooth(x, y)
         
        pot_analytic = potential_analytic(xx, yy)
        plt.plot(y, pot_analytic.T[nx/2 + nx/20, :], label='Analytic')
          
        pot_numeric = potential_function(xx, yy)
        plt.plot(y, pot_numeric.T[nx/2 + nx/20, :], label='Numeric')
        plt.plot(y, potential_function_smooth(0, y)[:, 0], label='Numeric_smooth')
        plt.legend(loc=0)
        plt.show()
         
        print np.sum(np.abs(pot_analytic.T[nx/2, :] - pot_numeric.T[nx/2, :]) * np.diff(x)[0])
        plt.plot(pot_analytic.T[nx/2 + nx/20., :] - pot_numeric.T[nx/2 + nx/20., :])
        plt.plot(pot_analytic.T[nx/2, :] - pot_numeric.T[nx/2, :])
        plt.show()        

#         plot.plot_planar_sensor(width=width, 
#                        pitch=width,
#                        thickness = thickness,
#                        n_pixel=n_pixel, 
#                        V_backplane=0, 
#                        V_readout=1,
#                        potential_function=potential_function,
#                        field_function=None,
#                        mesh=False)

if __name__ == "__main__":
    unittest.main()
    