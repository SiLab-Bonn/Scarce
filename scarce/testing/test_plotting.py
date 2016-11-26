import unittest
import os

from scarce import plot, fields


class TestPlotting(unittest.TestCase):

    def setUp(self):
        if os.getenv('TRAVIS', False):
            from xvfbwrapper import Xvfb
            self.vdisplay = Xvfb()
            self.vdisplay.start()

        # To have plt.show() non blocking
        import matplotlib.pyplot as p
        p.switch_backend('Agg')

    def tearDown(self):
        pass

    def test_plot_planar(self):
        ''' Check plotting of planar sensor.
        '''

        thickness = 200  # [um]
        width = 40  # [um]

        def potential_function(x, y):
            return fields.get_weighting_potential_analytic(x, y, D=thickness, S=width, is_planar=True)

        def field_function(x, y):
            return fields.get_weighting_field_analytic(x, y, D=thickness, S=width, is_planar=True)

        # Plot with analytical field function
        plot.plot_planar_sensor(potential_function=potential_function,
                                width=width,
                                pitch=width,
                                thickness=thickness,
                                n_pixel=1,
                                V_backplane=0,
                                V_readout=1,
                                field_function=field_function)

        # Plot without a field function
        plot.plot_planar_sensor(potential_function=potential_function,
                                width=width,
                                pitch=width,
                                thickness=thickness,
                                n_pixel=1,
                                V_backplane=0,
                                V_readout=1,
                                field_function=None)


if __name__ == "__main__":
    unittest.main()
