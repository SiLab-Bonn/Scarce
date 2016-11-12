import unittest
import numpy as np

from scarce.testing import tools
from scarce import fields


class Test(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_w_potential_analytic(self):
        ''' Check analytic weighting potential calculation
            for planar sensors.
        '''

        self.assertTrue(tools.check_with_fixture(fields.get_weighting_potential,
                                                 x=np.linspace(-200, 200, 200),
                                                 y=np.linspace(0, 250, 200),
                                                 D=[100, 150, 200, 250],
                                                 S=[5, 50, 100, 150, 250],
                                                 is_planar=True))

    def test_w_field_analytic(self):
        ''' Check analytic weighting field calculation
            for planar sensors.
        '''

        self.assertTrue(tools.check_with_fixture(fields.get_weighting_field,
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
        
                E_w_x, E_w_y = fields.get_weighting_field(xx, yy, D=thickness, S=width, is_planar=True)
                Phi_w = fields.get_weighting_potential(xx, yy, D=thickness, S=width, is_planar=True)
        
                E_w_x_2, E_w_y_2  = np.gradient(-Phi_w, np.gradient(x), np.gradient(y), axis=(1, 0))
                
                array_1 = np.array([E_w_x, E_w_y])
                array_2 = np.array([E_w_x_2, E_w_y_2])

                # Assert that less than 1 % of the field poinst have an error > 1%
                self.assertLess(tools.compare_arrays(array_1, array_2, threshold=0.01), 0.01) 



if __name__ == "__main__":
    unittest.main()
    