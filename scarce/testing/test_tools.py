import unittest
import os
import numpy as np

from scarce import tools, geometry, fields


class TestTools(unittest.TestCase):

    @classmethod
    def tearDownClass(cls):
        os.remove('tmp.sc')

    def test_save_and_load(self):
        ''' Check the saving and loading to disk.
        '''

        # Create data
        mesh = geometry.mesh_planar_sensor(
            n_pixel=9,
            width=50.,
            thickness=300.,
            resolution=50.,
            filename='planar_mesh_example.msh')
        potential = fields.calculate_planar_sensor_potential(mesh=mesh,
                                                             width=50.,
                                                             pitch=45.,
                                                             n_pixel=9,
                                                             thickness=300.,
                                                             n_eff=5e12,
                                                             V_bias=-160.,
                                                             V_readout=0.,
                                                             V_bi=1.5)
        min_x = float(mesh.getFaceCenters()[0, :].min())
        max_x = float(mesh.getFaceCenters()[0, :].max())
        description = fields.Description(potential,
                                         min_x=min_x,
                                         max_x=max_x,
                                         min_y=0,
                                         max_y=300.,
                                         nx=202, ny=200)

        # Force the creation of the potential and field functions
        description.get_field(0, 0)

        # Store and reload object
        tools.save(description, 'tmp.sc')
        description_2 = tools.load('tmp.sc')

        self.assertTrue(np.all(description.pot_data == description_2.pot_data))
        self.assertTrue(
            np.all(description.potential_grid == description_2.potential_grid))
        self.assertTrue(np.all(np.array(description.get_field(description._xx,
                                                              description._yy))
                               == np.array(description_2.get_field(description_2._xx,
                                                                   description_2._yy))))

if __name__ == "__main__":
    import logging
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    unittest.main()
