import unittest
import os

from scarce.examples import plot_properties
from scarce.examples import sensor_planar_weighting
from scarce.examples import sensor_3D_weighting
from scarce.examples import potential_1D
from scarce.examples import sensor_planar
from scarce.examples import sensor_3D
from scarce.examples import transient_planar
from scarce.examples import transient_3D


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

    def test_planar_sensor_weighting(self):
        ''' Check example to create planar weighting potential/field.
        '''
        sensor_planar_weighting.sensor_planar()

    def test_3D_sensor_weighting(self):
        ''' Check example to create 3D weighting potential/field.
        '''
        sensor_3D_weighting.sensor_3D()

    def test_1D_potential(self):
        ''' Check example to create a 1D potential of a planar sensor.
        '''
        potential_1D.create_1D_planar_sensor()

    def test_planar_sensor(self):
        ''' Check example to create a 2D planar sensor.
        '''
        sensor_planar.sensor_planar()

    def test_3D_sensor(self):
        ''' Check example to create 3D potential/field.
        '''
        sensor_3D.sensor_3D()

    def test_transient_planar(self):
        ''' Check example to calculate transient signal in planar sensor.
        '''
        transient_planar.transient_planar()

    def test_transient_3D(self):
        ''' Check example to calculate transient signal in 3D sensor.
        '''
        transient_3D.transient_3D()

if __name__ == "__main__":
    unittest.main()
