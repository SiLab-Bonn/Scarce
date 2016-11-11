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

#         tools.create_fixture(get_weighting_potential,
#                              x=np.linspace(-200, 200, 200),
#                              y=np.linspace(0, 250, 200),
#                              D=[100, 150, 200, 250],
#                              S=[5, 50, 100, 150, 250],
#                              is_planar=True)

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

#         tools.create_fixture(get_weighting_field,
#                              x=np.linspace(-200, 200, 200),
#                              y=np.linspace(0, 250, 200),
#                              D=[100, 150, 200, 250],
#                              S=[5, 50, 100, 150, 250],
#                              is_planar=True)

        self.assertTrue(tools.check_with_fixture(fields.get_weighting_field,
                                                 x=np.linspace(-200, 200, 200),
                                                 y=np.linspace(0, 250, 200),
                                                 D=[100, 150, 200, 250],
                                                 S=[5, 50, 100, 150, 250],
                                                 is_planar=True))

    @unittest.SkipTest
    def test_w_pot_field_analytic(self):
        ''' Check weighting potential/field of planar sensor.

        The check is if grad(-Phi) = E_w
        '''

        width = 100
        thickness = 200
        x = np.linspace(-width * 2, width * 2, 100)
        y = np.linspace(0, thickness, 100)
        xx, yy = np.meshgrid(x, y, sparse=True)

        e_w_x, e_w_y = fields.get_weighting_field(xx, yy, D=width, S=thickness, is_planar=True)
        p_w = fields.get_weighting_potential(xx, yy, D=width, S=thickness, is_planar=True)

        e_w_y_2, e_w_x_2 = np.gradient(-p_w, np.gradient(y), np.gradient(x))

        print e_w_x[10]
        print e_w_x_2[10]

        self.assertTrue(np.allclose(e_w_x, e_w_x_2))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    #     unittest.main()

    width = 100
    thickness = 200
    x = np.linspace(-width * 2, width * 2, 1000)
    y = np.linspace(0, thickness, 1000)
    xx, yy = np.meshgrid(x, y, sparse=True)

    e_w_x, e_w_y = fields.get_weighting_field(xx, yy, D=width, S=thickness, is_planar=True)
    p_w = fields.get_weighting_potential(xx, yy, D=width, S=thickness, is_planar=True)

    e_w_y_2, e_w_x_2 = np.gradient(-p_w, np.gradient(y), np.gradient(x))

    print e_w_x[10]
    print e_w_x_2[10]

#     self.assertTrue(np.allclose(e_w_x, e_w_x_2))
