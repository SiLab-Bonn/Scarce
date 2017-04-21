import unittest
import os

from scarce import plot, fields, geometry


class TestPlotting(unittest.TestCase):

    def setUp(self):
        if os.getenv('TRAVIS', False):
            from xvfbwrapper import Xvfb
            self.vdisplay = Xvfb()
            self.vdisplay.start()

        # To have plt.show() non blocking
        import matplotlib.pyplot as p
        p.switch_backend('Agg')

    @classmethod
    def tearDownClass(cls):
        os.remove('planar_mesh_tmp_3.msh')

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
        plot.plot_planar_sensor(pot_func=potential_function,
                                width=width,
                                pitch=width,
                                thickness=thickness,
                                n_pixel=1,
                                V_backplane=0,
                                V_readout=1,
                                field_func=field_function)

        # Plot without a field function
        plot.plot_planar_sensor(pot_func=potential_function,
                                width=width,
                                pitch=width,
                                thickness=thickness,
                                n_pixel=1,
                                V_backplane=0,
                                V_readout=1,
                                field_func=None)

    def test_plot_mesh(self):
        mesh = geometry.mesh_planar_sensor(
            n_pixel=5,
            width=50.,
            thickness=100.,
            resolution=100.,
            filename='planar_mesh_tmp_3.msh')

        plot.plot_mesh(mesh)

if __name__ == "__main__":
    import logging
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    unittest.main()
