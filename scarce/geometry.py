''' Function to geometrically describe and mesh planar and 3D electrode configurations. '''

import numpy as np
from scipy.interpolate import RectBivariateSpline, griddata

from fipy import GmshImporter2D
import pygmsh as pg
import meshio as mio


def mesh_3D_sensor(x, y, n_pixel_x, n_pixel_y, radius, nD, resolution):
    
    def origin_x_offsets(n_pixel_x, x):
        n_pixel_left = int((n_pixel_x) / 2)
        n_pixel_right = int((n_pixel_x - 1) / 2)
        
        for i in range(-n_pixel_left, n_pixel_right + 1, 1):
            yield i * x

    def generate_ro_pillars(geom, x, y, n_pixel_x, n_pixel_y, radius, nD, resolution, x0=0., y0=0.):
        pillars = []
        for offset_x in origin_x_offsets(n_pixel_x, x):
            pillars.extend(generate_ro_pillar_per_pixel(geom, x, y, radius, nD, resolution, x0 + offset_x, y0))
        return pillars

    def generate_ro_pillar_per_pixel(geom, x, y, radius, nD, resolution, x0, y0):
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
            [x0 - x / 2, y0 + y / 2, 0], lcar=resolution_x), points[2]]))
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
    
    def generate_edges_new(xl, xr, yb, yt, r):
        points = []
        # Left edge
        points.append(geom.add_point(
            [xl, yb + r, 0], lcar=resolution_x))
        points.append(geom.add_point(
            [xl, yt - r, 0], lcar=resolution_x))

        # Left, top
        points.append(geom.add_point(
            [xl + r, yt, 0], lcar=resolution_x))
        points.append(geom.add_point(
            [xr - r, yt, 0], lcar=resolution_x))

#         # Right top
#         points.append(
#             geom.add_point([xr, yt - r, 0], lcar=resolution_x))
#         points.append(
#             geom.add_point([xr, yt + r, 0], lcar=resolution_x))
# 
#         # Right edge
#         points.append(
#             geom.add_point([x0 + pitch_x / 2, y0 + pitch_y / 2 - r, 0], lcar=resolution_x))
#         points.append(
#             geom.add_point([x0 + pitch_x / 2, y0 + r - pitch_y / 2, 0], lcar=resolution_x))
# 
#         # Right bottom
#         points.append(
#             geom.add_point([x0 + pitch_x / 2 - r, y0 - pitch_y / 2, 0], lcar=resolution_x))
#         points.append(
#             geom.add_point([x0 + r, y0 - pitch_y / 2, 0], lcar=resolution_x))
# 
#         # Left bottom
#         points.append(
#             geom.add_point([x0 - r, y0 - pitch_y / 2, 0], lcar=resolution_x))
#         points.append(geom.add_point(
#             [x0 - (n_pixel_left + 1. / 2.) * pitch_x + r, y0 - pitch_y / 2, 0], lcar=resolution_x))

        return points


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
        
        n_pixel_left = int((n_pixel_x) / 2)
        n_pixel_right = int((n_pixel_x - 1) / 2)

        points = generate_edges(x, y,
                                n_pixel_x, n_pixel_y,
                                r, x0, y0)

        print 'points', points
        edge_pillars = generate_edge_pillars(points,
                                             x, y,
                                             n_pixel_x, n_pixel_y,
                                             x0, y0)
        pillars = generate_ro_pillars(geom,
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
