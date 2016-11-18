import logging
import fipy
import numpy as np
import meshio as mio

from scipy.interpolate import RectBivariateSpline, griddata

from scarce import silicon
from scarce import geometry
from scarce import plot


class Description(object):

    ''' Class to describe potential and field at any
        point in space. The numerical potential estimation
        is used and interpolated to make this possible. The field
        is derived from a smoothed potential interpolation to minimize
        numerical instabilities.
    '''

    def __init__(self, potential, min_x, max_x, min_y, max_y, nx=200, ny=200, smoothing=0.1):
        self.potential_data = potential
        self.potential_grid = self.interpolate_potential(self.potential_data)
        self.smoothing = smoothing

        # Do not calculate field on init, since it is time consuming
        # and maybe not needed
        self.potential_smooth = None
        self.field_x = None
        self.field_y = None

        self._x = np.linspace(min_x, max_x, nx)
        self._y = np.linspace(min_y, max_y, ny)

        # Create x,y plot grid
        self._xx, self._yy = np.meshgrid(self._x, self._y, sparse=True)

    def interpolate_potential(self, potential=None):
        ''' Interpolates the potential on a grid.
        '''
        if potential is None:
            potential = self.potential_data

        points = np.array(potential.mesh.getFaceCenters()).T
        values = np.array(potential.arithmeticFaceValue())

        def grid_interpolator(grid_x, grid_y):
            return griddata(points=points,
                            values=values,
                            xi=(grid_x, grid_y),
                            method='cubic',
                            rescale=False,
                            fill_value=np.nan)

        return grid_interpolator

    def get_potential(self, x, y):
        return self.potential_grid(x, y)

    def get_potential_smooth(self, x, y):
        if self.potential_smooth is None:
            self._smooth_potential()
        return self.potential_smooth(x, y).T

    def get_field(self, x, y):
        if self.field_x is None or self.field_y is None:
            self._derive_field()
        return np.array([self.field_x(x, y).T, self.field_y(x, y).T])

    def _smooth_potential(self, smoothing=None):
        ''' This function takes the potential grid interpolation
            and smooths the data points.

            Smoothing is really buggy in scipy, the only
            working way to smooth is apperently to smooth
            on a grid, thus mesh points of the potential
            solution cannot be used directly.
        '''

        if not smoothing:
            smoothing = self.smoothing

        def interpolate_nan(a):
            ''' Fills nans with closest non nan value.
            Might not work for multi dimensional arrays. :TODO:
            '''
            mask = np.isnan(a)
            a[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), a[~mask])
            return a

        # Interpolate potential on the grid
        potential_grid = self.potential_grid(self._xx, self._yy)

        # Fill nans otherwise Spline interpolations fails without error...
        potential_grid = interpolate_nan(potential_grid)

        # Smooth on the interpolated grid
        self.potential_smooth = RectBivariateSpline(self._xx, self._yy, potential_grid.T, s=smoothing, kx=3, ky=3)

    def _derive_field(self):
        ''' Takes the potential to calculate the field in x, y
        via E_x, E_y = - grad(Potential)
        with spline interpolation and smoothing.
        '''

        if not self.potential_smooth:
            self._smooth_potential()

        E_x, E_y = np.gradient(-self.potential_smooth(self._xx, self._yy), np.diff(self._x)[0], np.diff(self._y)[0])

        # Create spline interpolators for E_x,E_y
        self.field_x = RectBivariateSpline(self._x, self._y, E_x, s=0, kx=2, ky=2)
        self.field_y = RectBivariateSpline(self._x, self._y, E_y, s=0, kx=2, ky=2)


def calculate_planar_sensor_w_potential(mesh, width, pitch, n_pixel, thickness):
    logging.info('Calculating weighting potential')
    # Mesh validity check
    mesh_width = mesh.getFaceCenters()[0, :].max() - mesh.getFaceCenters()[0, :].min()

    if mesh_width != width * n_pixel:
        raise ValueError('The provided mesh width does not correspond to the sensor width')

    if mesh.getFaceCenters()[1, :].min() != 0:
        raise ValueError('The provided mesh does not start at 0.')

    if mesh.getFaceCenters()[1, :].max() != thickness:
        raise ValueError('The provided mesh does not end at sensor thickness.')

    potential = fipy.CellVariable(mesh=mesh, name='potential', value=0.)
    permittivity = 1.
    potential.equation = (fipy.DiffusionTerm(coeff=permittivity) == 0.)

    # Calculate boundaries
    backplane = mesh.getFacesTop()
    readout_plane = mesh.getFacesBottom()

    electrodes = readout_plane
    bcs = [fipy.FixedValue(value=0., faces=backplane)]
    X, _ = mesh.getFaceCenters()
    for pixel in range(n_pixel):
        pixel_position = width * (pixel + 1. / 2.) - width * n_pixel / 2.
        bcs.append(fipy.FixedValue(value=1. if pixel_position == 0. else 0.,
                                   faces=electrodes &
                                   (X > pixel_position - pitch / 2.) &
                                   (X < pixel_position + pitch / 2.)))

    potential.equation.solve(var=potential, boundaryConditions=bcs)
    return potential


def calculate_planar_sensor_potential(mesh, width, pitch, n_pixel, thickness,
                                      n_eff, V_backplane, V_readout=0):

    # Mesh validity check
    mesh_width = mesh.getFaceCenters()[0, :].max() - mesh.getFaceCenters()[0, :].min()
    if mesh_width != width:
        raise ValueError('The provided mesh width does not correspond to the sensor width (%d != %d)',
                         mesh_width,
                         width)

    if mesh.getFaceCenters()[1, :].min() != 0:
        raise ValueError('The provided mesh does not start at 0.')

    if mesh.getFaceCenters()[1, :].max() != thickness:
        raise ValueError('The provided mesh does not end at sensor thickness.')

    potential = fipy.CellVariable(mesh=mesh, name='potential', value=0.)
    permittivity = 1.
    potential.equation = (fipy.DiffusionTerm(coeff=permittivity) == 0.)

    # Calculate boundaries
    V_backplane = V_backplane
    backplane = mesh.getFacesTop()

    V_readout = V_readout
    readout_plane = mesh.getFacesBottom()

    electrodes = readout_plane
    bcs = [fipy.FixedValue(value=V_backplane, faces=backplane)]
    X, _ = mesh.getFaceCenters()
    for pixel in range(n_pixel):
        pixel_position = width * (pixel + 1. / 2.) - width * n_pixel / 2.
        bcs.append(fipy.FixedValue(value=V_readout,
                                   faces=electrodes &
                                   (X > pixel_position - pitch / 2.) &
                                   (X < pixel_position + pitch / 2.)))

    potential.equation.solve(var=potential, boundaryConditions=bcs)
    return potential


def get_weighting_potential_analytic(x, y, D, S, is_planar=True):
    """ Planar sensor:
        From Nuclear Instruments and Methods in Physics Research A 535 (2004)
        554-557, with correction from wbar = pi*w/2/D to wbar = pi*w/D with:

        x [um] is the offset from the middle of the electrode
        y [um] the position in the sensor
        D [um] the sensor thickness
        S [um] the pixel pitch

        3D sensor:
        Weighting potential for two cylinders with:
        D [um] distance between columns
        S [um] is the radius
    """

    # Wheighting potential for one pixel
    if is_planar:
        xbar = np.pi * x / D
        ybar = np.pi * (y - D) / D
        wbar = np.pi * S / D
        return -1. / np.pi * (np.arctan(np.tan(ybar / 2) * np.tanh((xbar + wbar / 2.) / 2.)) -
                              np.arctan(np.tan(ybar / 2) * np.tanh((xbar - wbar / 2.) / 2.)))
    else:
        R = S
        D = D / 2.  # D is the total distance between the columns
        a = np.sqrt(D * D - R * R)
        Phi_w = 1. / (4 * np.arccosh(D / R)) * \
            np.log(((x - a) ** 2 + y ** 2) / ((x + a) ** 2 + y ** 2)) + 0.5

        # Stability
        Phi_w = np.ma.masked_where(
            np.sqrt((x + D) * (x + D) + y * y) < R, Phi_w)
        Phi_w = np.ma.masked_where(
            np.sqrt((x - D) * (x - D) + y * y) < R, Phi_w)
        Phi_w = np.ma.masked_where(Phi_w < 0., Phi_w)
        Phi_w = np.ma.masked_where(Phi_w > 1., Phi_w)

        return Phi_w


def get_weighting_field_analytic(x, y, D, S, is_planar=True):
    """ From Nuclear Instruments and Methods in Physics Research A 535 (2004)
        554-557, with correction from wbar = pi*w/2/D to wbar = pi*w/D
        with x [um] is the position in the sensor [0:thickness], y [um] the offset from the
        middle of the electrode, D [um] the sensor thickness and S [um] the
        eletrode width. The field is calculated from the drivation of the
        potential in x and y.
    """

    if is_planar:
        xbar = np.pi * x / D
        ybar = np.pi * (y - D) / D
        wbar = np.pi * S / D

        # Not easy to find a more simple form
        denom = (np.cosh(1. / 2. * (wbar - 2. * xbar)) + np.cos(ybar)) * \
            (np.cosh(1. / 2. * (wbar + 2. * xbar)) + np.cos(ybar)) * D

        E_x = - np.sin(ybar) * np.sinh(wbar / 2.) * np.sinh(xbar) / denom

        E_y = np.sinh(
            wbar / 2.) * (np.cosh(wbar / 2.) + np.cos(ybar) * np.cosh(xbar)) / denom

        return E_x, E_y
    else:
        # 3D sensor:
        # From the analytical derivation of the get_weighting_potential function
        # Weighting potential for two cylinders with:
        # S [um] is the radius
        # D [um] distance between columns

        R = S
        D = D / 2.
        a = np.sqrt(D * D - R * R)

        E_x = a / (np.arccosh(D / R)) * (a ** 2 - x ** 2 + y ** 2) / \
            (((a - x) ** 2 + y ** 2) * ((a + x) ** 2 + y ** 2))
        E_y = -2 * a / (np.arccosh(D / R)) * (x * y) / \
            (((a - x) ** 2 + y ** 2) * ((a + x) ** 2 + y ** 2))

        E_x = np.ma.masked_where(np.sqrt((x + D) * (x + D) + y * y) < R, E_x)
        E_x = np.ma.masked_where(np.sqrt((x - D) * (x - D) + y * y) < R, E_x)
        E_y = np.ma.masked_where(np.sqrt((x + D) * (x + D) + y * y) < R, E_y)
        E_y = np.ma.masked_where(np.sqrt((x - D) * (x - D) + y * y) < R, E_y)

        return -E_x, -E_y


def get_potential_planar_analytic(x, V_bias, n_eff, D):
    """ Calculates the potential [V] in a planar sensor as a
        function of the position x between the electrodes [um],
        the bias voltage V_bias [V], the effective doping
        concentration n_eff [cm^-3] and the sensor Width D [um].
        The analytical function from the detector book p. 93 is used.
    """

    V_dep = silicon.get_depletion_voltage(n_eff, D)  # Depletion voltage

    a = (V_bias - V_dep) / D
    b = -2. * V_dep / (D ** 2)
    return (a - b / 2 * x) * x


def get_electric_field_analytic(x, y, V_bias, n_eff, D, S=None, is_planar=True):
    """ Calculates the 2D electric field E_x, E_y [V/um]

    Planar sensor:
        Calculates the field E_y[V/um], E_x = 0 in a planar sensor as a
        function of the position x between the electrodes [um],
        the bias voltage V_bias [V], the effective doping
        concentration n_eff [cm^-3] and the sensor Width D [um].
        The analytical function from the detector book p. 93 is used.

    3D sensor:
        Calculates the field E_x/E_y [V/um] in a 3d sensor as a function of the position
        x,y between the electrodes [um], the bias Voltage V_bias [V], the effective
        doping concentration n_eff [cm^-3], the electrode distance D [um] and radius R [um].
        So far the same field like the weighting field is used --> space charge is ignored.
    """

    if is_planar:
        if S:
            raise NotImplementedError(
                'The electrode width cannot be set, only full fill factor supported!')
        V_dep = silicon.get_depletion_voltage(n_eff, D)  # Depletion voltage
        a = (V_bias - V_dep) / D
        b = -2. * V_dep / (D ** 2)
        E_y = a - b * y
        return np.zeros_like(E_y), E_y
    else:
        E_x, E_y = get_weighting_field_analytic(x, y, D, S, is_planar=False)
        E_x = E_x * V_bias
        E_y = E_y * V_bias

        return E_x, E_y


def calculate_3D_sensor_potential(pitch_x, pitch_y, n_pixel_x, n_pixel_y, radius, resolution, V_readout, V_bias, nD=2):
    points, cells = geometry.mesh_3D_sensor(x=pitch_x,
                                            y=pitch_y,
                                            n_pixel_x=n_pixel_x,
                                            n_pixel_y=n_pixel_y,
                                            radius=radius,
                                            nD=nD,
                                            resolution=resolution)

    mio.write('sensor.msh', points, cells)
    mesh = fipy.GmshImporter2D('sensor.msh')

    plot.plot_mesh(mesh)

#     potential = fipy.CellVariable(mesh=mesh, name='potential', value=0.)
#     permittivity = 1.
#     potential.equation = (fipy.DiffusionTerm(coeff=permittivity) == 0.)
#
#     bcs = []
#     allfaces = mesh.getExteriorFaces()
#     X,Y =  mesh.getFaceCenters()
#
# Readout pillars
#     for pillar in range(nD):
#         position = pitch_x / nD * (pillar + 1. / 2.) - pitch_x / 2.
#         ring = allfaces & ( (X-position)**2+(Y)**2 < (radius)**2)
#         bcs.append(fipy.FixedValue(value=V_readout,faces=ring))
#
# Bias pillars
# Edges
#     positions = [(- pitch_x / 2., - pitch_y / 2.),
#                  (+ pitch_x / 2., - pitch_y / 2.),
#                  (+ pitch_x / 2., + pitch_y / 2.),
#                  (- pitch_x / 2., + pitch_y / 2.)]
# Sides
#     positions += [(0, - pitch_y / 2.),
#                  (0, + pitch_y / 2.)]
#
#     for pos_x, pos_y in positions:
#         ring = allfaces & ( (X-pos_x)**2+(Y-pos_y)**2 < (radius)**2)
#         bcs.append(fipy.FixedValue(value=V_bias, faces=ring))

# Calculate boundaries
#     p_pillars = mesh.getFaces()
#     n_pillars = mesh.getFacesTop()
#
#     electrodes = readout_plane
#     bcs = [fipy.FixedValue(value=V_backplane, faces=backplane)]
#
#     for pixel in range(n_pixel):
#         pixel_position = width * (pixel + 1. / 2.) - width * n_pixel / 2.
#         bcs.append(fipy.FixedValue(value=V_readout,
#                                    faces=electrodes &
#                                    (X > pixel_position - pitch / 2.) &
#                                    (X < pixel_position + pitch / 2.)))

#     potential.equation.solve(var=potential, boundaryConditions=bcs)
#     return potential

if __name__ == '__main__':
    pitch_x = 250.
    pitch_y = 50.
    n_pixel_x, n_pixel_y = 1, 1
    radius = 6.
    resolution = 50.
    V_readout, V_bias, = 0, -1

    potential = calculate_3D_sensor_potential(pitch_x, pitch_y, n_pixel_x, n_pixel_y, radius, resolution, V_readout, V_bias)
#     plot.plot_mesh(potential.mesh)
#     viewer = fipy.viewers.Viewer(vars=(potential, ))
#     viewer.plot("3D.png")

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
