''' Function to geometrically describe and mesh planar and 3D electrode configurations. '''

import numpy as np
from scipy.interpolate import SmoothBivariateSpline, interp2d
from scipy.interpolate import griddata

from fipy import GmshImporter2D
import pygmsh as pg
import meshio as mio
import tempfile


def mesh_3D_sensor(x, y, n_pixel_x, n_pixel_y, radius, nD, resolution):

    def generate_ro_pillar(geom, x, y, n_pixel_x, n_pixel_y, radius, nD, resolution, x0=0., y0=0.):
        pillars = []

        # Create readout pillars
        for pillar in range(nD):
            position = x / nD * (pillar + 1. / 2.) - x / 2.
            circle = geom.add_circle(x0=[position + x0, y0, 0.0],
                                     radius=radius,
                                     lcar=resolution / 4.,
                                     num_sections=4,
                                     # If compound==False, the section borders have to be points of the
                                     # discretization. If using a compound circle, they don't; gmsh can
                                     # choose by itself where to point the
                                     # circle points.
                                     compound=False
                                     )
            pillars.append(geom.add_line_loop(circle))

        return pillars

    def generate_edge_pillars(points, x, y, n_pixel_x, n_pixel_y, x0, y0):
        loop = []
        loop.append(geom.add_line(points[0], points[1]))
        loop.append(geom.add_circle_sector([points[1], geom.add_point(
            [x0 - x / 2, y0 + y / 2, 0], lcar=resolution_x / 4.), points[2]]))
        loop.append(geom.add_line(points[2], points[3]))
        loop.append(geom.add_circle_sector(
            [points[3], geom.add_point([x0, y0 + y / 2, 0], lcar=resolution_x), points[4]]))

        loop.append(geom.add_line(points[4], points[5]))
        loop.append(geom.add_circle_sector([points[5], geom.add_point(
            [x0 + x / 2, y0 + y / 2, 0], lcar=resolution_x), points[6]]))
        loop.append(geom.add_line(points[6], points[7]))
        loop.append(geom.add_circle_sector([points[7], geom.add_point(
            [x0 + x / 2, y0 - y / 2, 0], lcar=resolution_x), points[8]]))

        loop.append(geom.add_line(points[8], points[9]))
        loop.append(geom.add_circle_sector(
            [points[9], geom.add_point([x0, y0 - y / 2, 0], lcar=resolution_x), points[10]]))

        loop.append(geom.add_line(points[10], points[11]))
        loop.append(geom.add_circle_sector([points[11], geom.add_point(
            [x0 - x / 2, y0 - y / 2, 0], lcar=resolution_x), points[0]]))

        return geom.add_line_loop(loop)

    def generate_edges(pitch_x, pitch_y, n_pixel_x, n_pixel_y, r, x0, y0):
        points = []

        n_pixel_left = int((n_pixel_x) / 2)
        n_pixel_right = int((n_pixel_x - 1) / 2)

        print 'n_pixel_left', n_pixel_left
        print 'n_pixel_right', n_pixel_right

        # Left edge
        points.append(geom.add_point(
            [x0 - (n_pixel_left + 1. / 2.) * pitch_x, y0 + r - pitch_y / 2, 0], lcar=resolution_x))
        points.append(geom.add_point(
            [x0 - (n_pixel_left + 1. / 2.) * pitch_x, y0 + pitch_y / 2 - r, 0], lcar=resolution_x))

        # Left, top
        points.append(geom.add_point(
            [x0 + r - (n_pixel_left + 1. / 2.) * pitch_x, y0 + pitch_y / 2, 0], lcar=resolution_x))
        points.append(
            geom.add_point([x0 - r, y0 + pitch_y / 2, 0], lcar=resolution_x))

        # Right top
        points.append(
            geom.add_point([x0 + r, y0 + pitch_y / 2, 0], lcar=resolution_x))
        points.append(
            geom.add_point([x0 + pitch_x / 2 - r, y0 + pitch_y / 2, 0], lcar=resolution_x))

        # Right edge
        points.append(
            geom.add_point([x0 + pitch_x / 2, y0 + pitch_y / 2 - r, 0], lcar=resolution_x))
        points.append(
            geom.add_point([x0 + pitch_x / 2, y0 + r - pitch_y / 2, 0], lcar=resolution_x))

        # Right bottom
        points.append(
            geom.add_point([x0 + pitch_x / 2 - r, y0 - pitch_y / 2, 0], lcar=resolution_x))
        points.append(
            geom.add_point([x0 + r, y0 - pitch_y / 2, 0], lcar=resolution_x))

        # Left bottom
        points.append(
            geom.add_point([x0 - r, y0 - pitch_y / 2, 0], lcar=resolution_x))
        points.append(geom.add_point(
            [x0 - (n_pixel_left + 1. / 2.) * pitch_x + r, y0 - pitch_y / 2, 0], lcar=resolution_x))

        return points

    def generate_3D_pixel(geom, x, y, n_pixel_x, n_pixel_y, r, nD, resolution, x0=0., y0=0.):

        points = generate_edges(x, y,
                                n_pixel_x, n_pixel_y,
                                r, x0, y0)

        print 'points', points
        edge_pillars = generate_edge_pillars(points,
                                             x, y,
                                             n_pixel_x, n_pixel_y,
                                             x0, y0)
        pillars = generate_ro_pillar(geom,
                                     x, y,
                                     n_pixel_x, n_pixel_y,
                                     radius=r, nD=nD,
                                     resolution=resolution_x,
                                     x0=x0, y0=y0)

        geom.add_plane_surface([edge_pillars] + pillars)

        raw_codes = ['lc = %f;' % (resolution_x / 8.),
                     'Field[1] = Attractor;',
                     'Field[1].EdgesList = {c1, c2, c3, c4, c5, c6};'
                     'Field[1].NNodesByEdge = %d;' % resolution,
                     'Field[2] = MathEval;',
                     'Field[2].F = Sprintf(\"F1^3 + %g\", lc);',
                     'Background Field = 2;\n']
        geom.add_raw_code(raw_codes)

    if n_pixel_x < 1 or n_pixel_y < 1:
        raise RuntimeError(
            'Invalid parameter n_pixel_x, n_pixel_y = %d, %d' % (n_pixel_x, n_pixel_y))

    geom = pg.Geometry()
    resolution_x = x / resolution

    generate_3D_pixel(
        geom, x, y, n_pixel_x, n_pixel_y, radius, nD, resolution, x0=0, y0=0)

    return pg.generate_mesh(geom)


def mesh_planar_sensor(n_pixel, width, thickness, resolution=1., filename='sensor.msh'):  # TODO: size independent resolution parameter
    if n_pixel < 3:
        raise NotImplementedError('Less than 3 pixels result in wrong boundaries. Choose more pixels!')
    if not n_pixel % 2:
        raise NotImplementedError('Chose an odd pixel number (symmetry reasons)')

    geom = pg.Geometry()

    x = n_pixel * width

    resolution_x = x / resolution

    points_xyz = [
        [x / 2, thickness, 0],
        [x / 2, 0, 0],
        [width * 1.5, 0, 0],  # center 3 pixel region, top
        [-width * 1.5, 0, 0],  # center 3 pixel region, top
        [-x / 2, 0, 0],
        [-x / 2, thickness, 0],
        [-width * 1.5, thickness, 0],  # center 3 pixel region, bottom
        [width * 1.5, thickness, 0],  # center 3 pixel region, bottom
    ]

    # Decrease resolution to 1/4 for areas next to the 3 center pixel
    points = []
    for i, point_xyz in enumerate(points_xyz):
        points.append(geom.add_point(point_xyz, lcar=resolution_x if i in (2, 3, 6, 7) else resolution_x * 4.))

    # Create lines
    lines = [geom.add_line(points[i], points[i + 1])
             for i in range(len(points) - 1)]
    lines.append(geom.add_line(points[-1], points[0]))

    line_loop = geom.add_line_loop(lines)
    geom.add_plane_surface([line_loop])

    # Add 1/x1.5 law for the mesh size
    raw_codes = ['lc = %f;' % (resolution_x / 4.),
                 'Field[1] = Attractor;',
                 'Field[1].EdgesList = {l3};'
                 'Field[1].NNodesByEdge = %d;' % resolution,
                 'Field[2] = MathEval;',
                 'Field[2].F = Sprintf(\"F1^3 + %g\", lc);',
                 'Background Field = 2;\n']

    geom.add_raw_code(raw_codes)

    points, cells = pg.generate_mesh(geom)

    mio.write(filename, points, cells)
    return GmshImporter2D(filename)


def interpolate_potential_old_smooth(potential, smoothing=None):
    x = np.array(potential.mesh.getFaceCenters()[0])
    y = np.array(potential.mesh.getFaceCenters()[1])
    z = np.array(potential.arithmeticFaceValue())
    return SmoothBivariateSpline(x, y, z, s=smoothing, kx=3, ky=3)


def interpolate_potential_old_smooth_2(potential):
    x = np.array(potential.mesh.getFaceCenters()[0])
    y = np.array(potential.mesh.getFaceCenters()[1])
    z = np.array(potential.arithmeticFaceValue())
    return interp2d(x, y, z, kind='cubic')


def interpolate_potential(potential):
    points = np.array(potential.mesh.getFaceCenters()).T
    values = np.array(potential.arithmeticFaceValue())

    def interpolator(grid_x, grid_y):
        return griddata(points=points,
                        values=values,
                        xi=(grid_x, grid_y),
                        method='cubic',
                        rescale=False)

    return interpolator


def calculate_field(potential):
    ''' Takes the potential to calculate the field in x, y
    with slight smoothing via: E_x, E_y = - grad(Potential)
    '''

    field = potential.faceGradAverage

    mesh = field.mesh

    C = mesh.faceCenters

    def interpolator(grid_x, grid_y):
        U = griddata(C.value.T, field.value[0],
                     (grid_x, grid_y), method='cubic')
        V = griddata(C.value.T, field.value[1],
                     (grid_x, grid_y), method='cubic')
        return -U, -V

    return interpolator


if __name__ == '__main__':
    from scarce import fields, plot
#     pitch_x = 250.
#     pitch_y = 50.
#     n_pixel_x, n_pixel_y = 1, 1
#     radius = 6.
#     resolution = 50.
#     V_readout, V_bias,  = 0, -1
#
#     potential = calculate_3D_sensor_potential(pitch_x, pitch_y, n_pixel_x, n_pixel_y, radius, resolution, V_readout, V_bias)
# #     plot.plot_mesh(potential.mesh)
# #     viewer = fipy.viewers.Viewer(vars=(potential, ))
# #     viewer.plot("3D.png")
#
#     min_x, max_x = np.min(np.array(potential.mesh.getFaceCenters()[0])), np.max(np.array(potential.mesh.getFaceCenters()[0]))
#     min_y, max_y = np.min(np.array(potential.mesh.getFaceCenters()[1])), np.max(np.array(potential.mesh.getFaceCenters()[1]))
#
#     print 'Interpolate'
#
#     xnew = np.linspace(min_x, max_x, 1000)
#     ynew = np.linspace(min_y, max_y, 1000)
#     xnew_plot, ynew_plot = np.meshgrid(xnew, ynew)
#
#     potential_function = interpolate_potential_2(potential)
#     print 'Done'
#
#     plot.plot_3D_sensor(potential_function,
#                         pitch_x,
#                         pitch_y,
#                         n_pixel,
#                         radius,
#                         V_bias,
#                         V_readout,
#                         min_x,
#                         max_x,
#                         min_y,
#                         max_y
#                         )

    width = 200
    pitch = 240
    n_pixel = 1
    thickness = 250
    resolution = 50
    V_backplane, V_readout = -1, 0

    mesh = mesh_planar_sensor(x=width * n_pixel,
                              thickness=thickness,
                              resolution=resolution)
    potential = fields.calculate_planar_sensor_potential(mesh, width, pitch, n_pixel, thickness, V_backplane, V_readout)

    plot.plot_mesh(mesh, invert_y_axis=True)

#     plot.plot_mesh(potential.mesh, invert_y_axis=True)

#     min_x, max_x = np.min(np.array(potential.mesh.getFaceCenters()[0])), np.max(np.array(potential.mesh.getFaceCenters()[0]))
#     min_y, max_y = np.min(np.array(potential.mesh.getFaceCenters()[1])), np.max(np.array(potential.mesh.getFaceCenters()[1]))
#
#     print 'Interpolate', np.square(abs(V_backplane - V_readout))
#     potential_function = interpolate_potential(potential)
#     plot.plot_planar_sensor(potential_function,
#                             width,
#                             pitch,
#                             n_pixel,
#                             thickness,
#                             V_backplane,
#                             V_readout,
#                             min_x,
#                             max_x,
#                             min_y,
#                             max_y)
