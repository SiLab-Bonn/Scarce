''' gmsh has to be installed and put into the ENVIRONMENT variable in a way that
    the command gmsh from the terminal/console starts it.
'''

import pygmsh as pg


def generate(resolution):
    geom = pg.Geometry()

    circle = geom.add_circle(x0=[0.0, 0.0, 0.0],
                             radius=1.0,
                             lcar=resolution,
                             num_sections=4,
                             # If compound==False, the section borders have to be points of the
                             # discretization. If using a compound circle, they don't; gmsh can
                             # choose by itself where to point the circle points.
                             compound=True
                             )

    ll = geom.add_line_loop(circle)
    geom.add_plane_surface(ll)

    return geom


if __name__ == '__main__':
    points, cells = pg.generate_mesh(generate(resolution=0.3))
    import matplotlib.pyplot as plt

    # Plot triangle cells
    plt.plot(points[cells['triangle'].T, 0], points[cells['triangle'].T, 1], '-', color='black')
    plt.xlim((-1.1, 1.1))
    plt.ylim((-1.1, 1.1))
    plt.show()

    # Save mesh
    # import meshio as mio
    # mio.write('circle.vtk', points, cells)
