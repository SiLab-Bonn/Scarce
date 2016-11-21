''' Function to geometrically describe and mesh planar and 3D electrode configurations. '''

import logging

from fipy import GmshImporter2D
import pygmsh as pg
import meshio as mio


class SensorDescription3D(object):

    """ This class provides helper functions to describe a 3D sensor array.

    """

    def __init__(self, width_x, width_y, n_pixel_x, n_pixel_y, radius, nD, x0=0., y0=0.):
        if n_pixel_x < 1 or n_pixel_y < 1:
            raise RuntimeError(
                'Invalid parameter n_pixel_x, n_pixel_y = %d, %d' % (n_pixel_x, n_pixel_y))

        self.width_x = width_x
        self.width_y = width_y
        self.n_pixel_x = n_pixel_x
        self.n_pixel_y = n_pixel_y
        self.radius = radius
        self.nD = nD
        self.x0 = x0
        self.y0 = y0

    def get_n_pixel_left(self):
        return int((self.n_pixel_x) / 2)

    def get_n_pixel_right(self):
        return int((self.n_pixel_x - 1) / 2)

    def get_n_pixel_top(self):
        return int((self.n_pixel_y) / 2)

    def get_n_pixel_bottom(self):
        return int((self.n_pixel_y - 1) / 2)

    def get_array_corners(self):
        min_x = self.x0 - self.width_x / 2. - self.get_n_pixel_left() * self.width_x
        max_x = self.x0 + self.width_x / 2. + self.get_n_pixel_right() * self.width_x
        min_y = self.y0 - self.width_y / 2. - self.get_n_pixel_top() * self.width_y
        max_y = self.y0 + self.width_y / 2. + self.get_n_pixel_bottom() * self.width_y

        return min_x, max_x, min_y, max_y

    def get_pixel_x_offsets(self):
        ''' Pixel center offsets in x for the number of pixels in x '''

        for i in range(-self.get_n_pixel_left(), self.get_n_pixel_right() + 1, 1):
            yield i * self.width_x

    def get_pixel_y_offsets(self):
        ''' Pixel center offsets in y for the number of pixels in y '''

        for i in range(-self.get_n_pixel_top(), self.get_n_pixel_bottom() + 1, 1):
            yield i * self.width_y

    def get_side_col_x_offsets(self):
        ''' Pixel side column offsets in x excluding edges '''

        offsets = []
        for ix0 in self.get_pixel_x_offsets():
            for pillar in range(self.nD):
                offsets.append(self.width_x / self.nD * (pillar + 1.) - self.width_x / 2. + ix0)
        return offsets[:-1]  # Last offset is edge column

    def get_side_col_y_offsets(self):
        ''' Pixel side column offsets in y excluding edges '''

        offsets = []
        for iy0 in self.get_pixel_y_offsets():
            offsets.append(self.width_y / 2. + iy0)
        return offsets[:-1]  # Last offset is edge column

    def get_ro_col_offsets(self):
        ''' Offsets of the readout columns in x,y '''

        for offset_x in self.get_pixel_x_offsets():
            for offset_y in self.get_pixel_y_offsets():
                for pillar in range(self.nD):
                    yield self.width_x / self.nD * (pillar + 1. / 2.) - self.width_x / 2. + offset_x, offset_y

    def get_center_bias_col_offsets(self):
        ''' Offsets of the center bias columns in x,y '''

        for offset_x in self.get_side_col_x_offsets():
            for offset_y in list(self.get_pixel_y_offsets())[:-1]:
                yield offset_x, offset_y + self.width_y / 2.

    def get_side_bias_col_offsets(self):
        ''' Offsets of the side bias columns in x,y '''

        offsets_x, offsets_y = [], []
        min_x, max_x, min_y, max_y = self.get_array_corners()

        # Left side
        for y in self.get_side_col_y_offsets():
            offsets_x.append(min_x)
            offsets_y.append(y)

        # Right side
        for y in self.get_side_col_y_offsets():
            offsets_x.append(max_x)
            offsets_y.append(y)

        # Bottom side
        for x in self.get_side_col_x_offsets():
            offsets_x.append(x)
            offsets_y.append(min_y)

        # Top side
        for x in self.get_side_col_x_offsets():
            offsets_x.append(x)
            offsets_y.append(max_y)

        return zip(offsets_x, offsets_y)

    def get_edge_bias_col_offsets(self):
        ''' Offsets of the side bias columns in x,y '''

        offsets_x, offsets_y = [], []
        min_x, max_x, min_y, max_y = self.get_array_corners()

        offsets_x.append(min_x)
        offsets_y.append(min_y)

        offsets_x.append(max_x)
        offsets_y.append(min_y)

        offsets_x.append(min_x)
        offsets_y.append(max_y)

        offsets_x.append(max_x)
        offsets_y.append(max_y)

        return zip(offsets_x, offsets_y)

    def position_in_center_pixel(self, x, y):
        return -self.width_x / 2. <= x + self.x0 <= self.width_x / 2. and -self.width_y / 2. <= y + self.y0 <= self.width_y / 2.


def mesh_3D_sensor(width_x, width_y, n_pixel_x, n_pixel_y, radius, nD, resolution, filename='sensor.msh'):
    ''' Create the mesh of a 3D sensor array '''

    logging.info('Mesh 3D sensor array')

    desc = SensorDescription3D(width_x, width_y, n_pixel_x, n_pixel_y, radius, nD)

    # Center pixel column segments that need higher resolution
    high_res_segments = []

    def generate_ro_pillars(geom, radius, resolution, x0=0., y0=0.):
        pillars = []

        for position_x, position_y in desc.get_ro_col_offsets():
            circle = geom.add_circle(x0=[position_x + x0, position_y + y0, 0.0],
                                     radius=radius,
                                     lcar=resolution,
                                     num_sections=4,
                                     # If compound==False, the section borders have to be points of the
                                     # discretization. If using a compound circle, they don't; gmsh can
                                     # choose by itself where to point the
                                     # circle points.
                                     compound=False
                                     )
            # Check if circle belongs to center pixel
            if desc.position_in_center_pixel(position_x, position_y):
                high_res_segments.extend(circle)

            pillars.append(geom.add_line_loop(circle))

        return pillars

    def generate_bias_pillars(geom, radius, resolution, x0=0., y0=0.):
        pillars = []

        for position_x, position_y in desc.get_center_bias_col_offsets():
            circle = geom.add_circle(x0=[position_x + x0, position_y + y0, 0.0],
                                     radius=radius,
                                     lcar=resolution,
                                     num_sections=4,
                                     # If compound==False, the section borders have to be points of the
                                     # discretization. If using a compound circle, they don't; gmsh can
                                     # choose by itself where to point the
                                     # circle points.
                                     compound=False
                                     )

            # Check if circle belongs to center pixel
            if desc.position_in_center_pixel(position_x, position_y):
                high_res_segments.extend(circle)

            pillars.append(geom.add_line_loop(circle))

        return pillars

    def generate_domain_edge(x0_1, x0_2, y0_1, y0_2, r):
        ''' Generates the perimeter of the pixel array including 1/4 edge and
        1/2 side columns.
        '''

        points = []
        loop = []

        def add_partial_column(x0, y0, last=False):
            circle = geom.add_circle_sector(
                [points[-1] if last else points[-2], geom.add_point([x0, y0, 0], lcar=resolution_x), points[0] if last else points[-1]])

            # Check if circle belongs to center pixel
            if desc.position_in_center_pixel(x0, y0):
                high_res_segments.append(circle)

            loop.append(circle)

        def add_side_column(x0_1, x0_2, y0_1, y0_2):
            ''' Adds a halve circle at the x side of the pixel array.
            x0_1, x0_2, y0_1, y0_2 are start/end points in x and y.
            '''

            points.append(geom.add_point(
                [x0_1, y0_1, 0], lcar=resolution_x))
            loop.append(geom.add_line(points[-2], points[-1]))
            points.append(
                geom.add_point([x0_2, y0_2, 0], lcar=resolution_x))
            add_partial_column(x0=(x0_1 + x0_2) / 2., y0=(y0_1 + y0_2) / 2.)

        # Left side
        points.append(geom.add_point(
            [x0_1, y0_1 + r, 0], lcar=resolution_x))

        # Left side column halves
        for iy0 in desc.get_side_col_y_offsets():
            add_side_column(x0_1=x0_1, x0_2=x0_1, y0_1=iy0 - r, y0_2=iy0 + r)

        points.append(geom.add_point(
            [x0_1, y0_2 - r, 0], lcar=resolution_x))
        loop.append(geom.add_line(points[-2], points[-1]))

        # Left edge, bottom
        points.append(geom.add_point(
            [x0_1 + r, y0_2, 0], lcar=resolution_x))
        add_partial_column(x0=x0_1, y0=y0_2)

        # Bottom side column halve(s)
        for ix0 in desc.get_side_col_x_offsets():
            add_side_column(x0_1=ix0 - r, x0_2=ix0 + r, y0_1=y0_2, y0_2=y0_2)

        # Close bottom side
        points.append(geom.add_point(
            [x0_2 - r, y0_2, 0], lcar=resolution_x))
        loop.append(geom.add_line(points[-2], points[-1]))

        # Right edge, bottom
        points.append(
            geom.add_point([x0_2, y0_2 - r, 0], lcar=resolution_x))
        add_partial_column(x0=x0_2, y0=y0_2)

        # Right side column halve(s)
        for iy0 in desc.get_side_col_y_offsets()[::-1]:
            add_side_column(x0_1=x0_2, x0_2=x0_2, y0_1=iy0 + r, y0_2=iy0 - r)

        # Close right side
        points.append(
            geom.add_point([x0_2, y0_1 + r, 0], lcar=resolution_x))
        loop.append(geom.add_line(points[-2], points[-1]))

        # Right edge, top
        points.append(
            geom.add_point([x0_2 - r, y0_1, 0], lcar=resolution_x))
        add_partial_column(x0=x0_2, y0=y0_1)

        # Top side column halve(s)
        for ix0 in desc.get_side_col_x_offsets()[::-1]:  # return order
            add_side_column(x0_1=ix0 + r, x0_2=ix0 - r, y0_1=y0_1, y0_2=y0_1)

        # Close top side
        points.append(geom.add_point(
            [x0_1 + r, y0_1, 0], lcar=resolution_x))
        loop.append(geom.add_line(points[-2], points[-1]))

        # Right edge, top; closes loop
        add_partial_column(x0=x0_1, y0=y0_1, last=True)

        return [geom.add_line_loop(loop)]

    def generate_3D_array(geom, width_x, width_y, r, resolution, x0=0., y0=0.):
        min_x, max_x, min_y, max_y = desc.get_array_corners()

        # Perimeter of 3D array, including edge and side columns
        perimeter = generate_domain_edge(x0_1=min_x,
                                         x0_2=max_x,
                                         y0_1=min_y,
                                         y0_2=max_y,
                                         r=r)

        ro_pillars = generate_ro_pillars(geom,
                                         radius=r,
                                         resolution=resolution_x,
                                         x0=x0, y0=y0)

        bias_pillars = generate_bias_pillars(geom,
                                             radius=r,
                                             resolution=resolution_x,
                                             x0=x0, y0=y0)

        geom.add_plane_surface(perimeter + ro_pillars + bias_pillars)

        raw_codes = ['lc = %f;' % (resolution_x / 8.),
                     'Field[1] = Attractor;',
                     'Field[1].EdgesList = {%s};' % ','.join(high_res_segments),
                     'Field[1].NNodesByEdge = %d;' % resolution,
                     'Field[2] = MathEval;',
                     'Field[2].F = Sprintf(\"F1^3 + %g\", lc);',
                     'Background Field = 2;\n']
        geom.add_raw_code(raw_codes)

    geom = pg.Geometry()
    resolution_x = width_x / resolution

    generate_3D_array(
        geom, width_x, width_y, radius, resolution, x0=0., y0=0.)

    points, cells = pg.generate_mesh(geom)

    mio.write(filename, points, cells)
    return GmshImporter2D(filename)


def mesh_planar_sensor(n_pixel, width, thickness, resolution=1., filename='sensor.msh'):  # TODO: size independent resolution parameter
    if n_pixel < 3:
        raise logging.warning('Less than 3 pixels result in quite wrong boundaries. It is better to choose more pixelss!')
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
