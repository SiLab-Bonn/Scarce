import unittest
import os
import meshio

from scarce import geometry


class TestGeometry(unittest.TestCase):

    @classmethod
    def tearDownClass(cls):
        os.remove('planar_mesh_tmp.msh')

    def test_mesh_planar(self):
        ''' Check the mesh generation for planar sensor.

        Likely due to random pertubation (on purpose) the result
        is never exactly constant. Thus only numbers of points checked.
        '''

        n_pixel = 5
        width = 50
        thickness = 200

        geometry.mesh_planar_sensor(n_pixel,
                                    width,
                                    thickness=thickness,
                                    resolution=100,
                                    filename='planar_mesh_tmp.msh')

        points, _, _, _, _ = meshio.read('planar_mesh_tmp.msh')

        self.assertGreater(points.shape[0], 14000)

    def test_mesh_3D(self):
        ''' Check for no excpetion during meshing. '''

        for n_pixel_x in [1, 2, 3]:
            for n_pixel_y in [1, 2, 3]:
                for nD in [1, 2, 3]:
                    geometry.mesh_3D_sensor(width_x=250,
                                            width_y=50,
                                            n_pixel_x=n_pixel_x,
                                            n_pixel_y=n_pixel_y,
                                            radius=6.,
                                            nD=nD,
                                            resolution=10)  # Low res for lower time

    def test_potential_interpolation(self):
        pass

if __name__ == "__main__":
    import logging
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    unittest.main()
