import unittest
import numpy as np

from siliconproperties.python_files.getWeightingField import get_weighting_field
from siliconproperties.python_files.getWeightingPotential import get_weighting_potential


class Test(unittest.TestCase):


    def setUp(self):
        pass


    def tearDown(self):
        pass


    def test_w_field_analytic(self):
        ''' Check weighting potential/field of planar sensor.
        
        The check is if grad(-Phi) = E_w
        '''

        width = 100
        thickness = 200
        x = np.linspace(-width * 2, width * 2, 100)
        y = np.linspace(0, thickness, 100)
        xx, yy = np.meshgrid(x, y, sparse=True)

        e_w = get_weighting_field(xx, yy, D=width, S=thickness, is_planar=True)
        p_w = get_weighting_potential(xx, yy, D=width, S=thickness, is_planar=True)
        
        self.assertTrue(np.allclose(e_w, np.gradient(-p_w)))
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()