import unittest
import os

from scarce.examples import plot_properties
from scarce.examples import sensor_planar_weighting
from scarce.examples import sensor_3D_weighting
from scarce.examples import potential_1D


class TestExamples(unittest.TestCase):

    def setUp(self):
        if os.getenv('TRAVIS', False):
            from xvfbwrapper import Xvfb
            self.vdisplay = Xvfb()
            self.vdisplay.start()

        # To have plt.show() non blocking
        import matplotlib.pyplot as p
        p.switch_backend('Agg')

    def tearDown(self):
        filelist = [f for f in os.listdir(".") if f.endswith(".pdf")]
        for f in filelist:
            os.remove(f)

    def test_plot_properties(self):
        ''' Check example to plot all define silicon properties.
        '''
        plot_properties.create_plots()

    def test_planar_sensor(self):
        ''' Check example to create planar weightinf potential/field.
        '''
        sensor_planar_weighting.sensor_planar()

    def test_3D_sensor(self):
        ''' Check example to create 3D weighting potential/field.
        '''
        sensor_3D_weighting.sensor_3D()

    def test_1D_potential(self):
        ''' Check example to create a 1D potential of a planar sensor.
        '''
        potential_1D.create_1D_planar_sensor()

if __name__ == "__main__":
    unittest.main()
