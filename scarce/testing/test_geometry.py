import unittest
import os
import meshio

from scarce import geometry


class Test(unittest.TestCase):

    @classmethod
    def tearDownClass(cls):
        os.remove('planar_mesh_tmp.msh')

    def test_mesh_planar(self):
        ''' Check the mesh generation for planar sensor.
        Likely due to random pertubation (on purpose) the result
        is never exactly constant.
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


if __name__ == "__main__":
    unittest.main()
