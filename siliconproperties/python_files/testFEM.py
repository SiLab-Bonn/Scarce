''' gmsh has to be installed and put into an environment variable (e.g. PATH)
    in a way that the command gmsh from the terminal/console starts it.
'''

import pygmsh as pg
import numpy as np


def generate(resolution):
    geom = pg.Geometry()

    # Two readout pillars
    circle_left = geom.add_circle(x0=[-1.5, 0.0, 0.0],
                                  radius=0.5,
                                  lcar=resolution / 2.,
                                  num_sections=4,
                                  # If compound==False, the section borders have to be points of the
                                  # discretization. If using a compound circle, they don't; gmsh can
                                  # choose by itself where to point the circle points.
                                  compound=False
                                  )
    circle_left_ll = geom.add_line_loop(circle_left)

    circle_right = geom.add_circle(x0=[1.5, 0.0, 0.0],
                                   radius=0.5,
                                   lcar=resolution / 2.,
                                   num_sections=4,
                                   # If compound==False, the section borders have to be points of the
                                   # discretization. If using a compound circle, they don't; gmsh can
                                   # choose by itself where to point the circle points.
                                   compound=False
                                   )
    circle_right_ll = geom.add_line_loop(circle_right)

    geom.add_rectangle(xmin=-3,
                       xmax=3,
                       ymin=-3,
                       ymax=3,
                       z=0,
                       lcar=resolution,
                       holes=[circle_left_ll, circle_right_ll])

    return geom


if __name__ == '__main__':
    points, cells = pg.generate_mesh(generate(resolution=0.3))
    import matplotlib.pyplot as plt

    # Plot triangle cells
    plt.plot(points[cells['triangle'].T, 0], points[cells['triangle'].T, 1], '-', color='black')
    plt.plot(points[cells['line'].T, 0], points[cells['line'].T, 1], '-', color='black')
#     plt.xlim((-2.1, 2.1))
#     plt.ylim((-2.1, 2.1))
    plt.show()

    # Save mesh
    # import meshio as mio
    # mio.write('circle.vtk', points, cells)
