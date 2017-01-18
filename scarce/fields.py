r"""The field module calculates the potentials and fields of silicon sensors.

The potential is determined numerically by solving these equations on a mesh:

.. math:: \nabla^2 \Phi = 0
   :label: laplace

.. math:: \nabla^2 \Phi = \frac{\rho}{\epsilon}
   :label: poisson

For the weighting potential equation :eq:`laplace` is solved with the
boundary conditions:

.. math::

   \begin{eqnarray}
      \Phi    & = & 0 \\
      \Phi_r    & = & 1
   \end{eqnarray}

The pixel readout electrode(s) are at a potential 1 and all other equipotential
pixel parts (backside, bias columns, etc.) at 0.

For the electric potential the equation :eq:`poisson` is solved with the
boundary conditions:

.. math::

   \begin{eqnarray}
      \Phi    & = & V_{bias} \\
      \Phi_r    & = & V_{readout}
   \end{eqnarray}

The pixel readout electrode(s) are at :math:`V_{readout}` potential and the bias parts
(backside, bias columns, etc.) are at :math:`V_{bias}`.

The field is then derived via:

.. math::
   \vec{E} = -\nabla \phi

.. NOTE::
   For simple cases (e.g. planar sensor with 100% fill factor) also analytical
   solutions are provided. The analytical results are also used to benchmark
   the numerical results in the automated unit tests.

"""

import fipy
import numpy as np
import logging

from scipy.interpolate import interp1d, RectBivariateSpline, griddata, interp2d, SmoothBivariateSpline
from scipy import constants

from scarce import silicon
from scarce import geometry
from scarce import constant as C
from scarce import solver

_LOGGER = logging.getLogger(__name__)


class Description(object):

    ''' Class to describe potential and field at any
        point in space. The numerical potential estimation
        is used and interpolated. The field is derived from
        a smoothed potential interpolation to minimize
        numerical instabilities.
    '''

    def __init__(self, potential, min_x, max_x, min_y, max_y,
                 nx, ny, smoothing=0.1):
        _LOGGER.debug('Create potential and field description')
        self.pot_data = potential

        try:
            self.depletion_data = np.array(
                [self.pot_data.depletion[0], self.pot_data.depletion[1]])
        except AttributeError:
            self.depletion_data = None

        self.min_x = min_x
        self.min_y = min_y
        self.max_x = max_x
        self.max_y = max_y
        self._nx = nx
        self._ny = ny

        self.potential_grid_inter = self.interpolate_potential(self.pot_data)
        self.smoothing = smoothing

        # Do not calculate field on init, since it is time consuming
        # and maybe not needed
        self.pot_smooth = None
        self.field_x = None
        self.field_y = None

        # Do not calculate depletion boundaries on init
        # since it is time consuming and maybe not needed
        self.depletion_region = None

        self._x = np.linspace(self.min_x, self.max_x, nx)
        self._y = np.linspace(self.min_y, self.max_y, ny)

        # Create sparse x,y plot grid
        self._xx, self._yy = np.meshgrid(self._x, self._y, sparse=True)

        # Potential interpolated on a grid with NaN set to closest value
        self.potential_grid = self.potential_grid_inter(self._xx, self._yy)
#         self._extrapolate_boundary()
        self.potential_grid = self._interpolate_nan(self.potential_grid)

    def _extrapolate_boundary(self):
        ''' Extrapolate the potental at the boundary with constant gradient
            to prevent oscillations when smoothing.
        '''

        # Extend array by 10% in y direction, y direction is
        # in dimension 0 here
        y_shape = int(self.potential_grid.shape[0] * 1.0)
        # Difference to old shape
        dy_shape = y_shape - self.potential_grid.shape[0]

        # Create extended array and fill with old data
        new_potential_grid = np.zeros(
            shape=(y_shape, self.potential_grid.shape[1]))  # ax2.plot(desc._y, -np.gradient(desc.potential_grid.T[xi, :], desc._y[1]-desc._y[0], edge_order=1), '-', label='DIFF N')
        new_potential_grid[dy_shape:, :] = self.potential_grid

        # Calculate the gradient in y for linear extrapolation
        dy = self.potential_grid[0, :] - self.potential_grid[1, :]

        # Extrapolate at the extended area
        slope = (np.tile(np.arange(1, dy_shape + 1)[::-1],
                         (self.potential_grid.shape[1], 1))).T * dy[np.newaxis, :]
        new_potential_grid[:dy_shape] = self.potential_grid[0, :] + slope
        # [:new_potential_grid.shape[0] - dy_shape, :]
        self.potential_grid = new_potential_grid

        # Also extend y positions
        self._y = np.linspace(
            self.min_y - dy_shape * np.diff(self._y)[0], self.max_y,
            self._ny + dy_shape)

        print 'self._y', self._y
#         raise

    def interpolate_potential(self, potential=None):
        ''' Interpolates the potential on a grid.
        '''
        _LOGGER.debug('Interpolate potential')
        if potential is None:
            potential = self.pot_data

        points = np.array(potential.mesh.getFaceCenters()).T
        values = np.array(np.array(potential.arithmeticFaceValue()))

        def grid_interpolator(grid_x, grid_y):
            return griddata(points=points,
                            values=values,
                            xi=(grid_x, grid_y),
                            method='linear',
                            rescale=False,
                            fill_value=np.nan)  # values.max())

        return grid_interpolator

    def get_potential_minimum(self, axis=None):
        ''' Returns the minimum potential value
        '''
        return self.potential_grid.min(axis=axis)

    def get_potential_minimum_pos_y(self):
        ''' Returns the position in the array with
            the potential minimum
        '''

        return self._y[np.argmin(self.potential_grid, axis=0)]

    def get_depletion(self, x):
        ''' Returns the depletion boundary at x.
            For planar sensors only!
        '''
        if not self.depletion_region:  # Calculate on demand to safe time
            _LOGGER.debug('Calculate depletion region description')
            if self.depletion_data is not None:
                self.depletion_region = interp1d(x=self.depletion_data[0],
                                                 y=self.depletion_data[1],
                                                 kind='cubic'
                                                 )
            else:
                raise RuntimeError(
                    'The data does not have depletion information.')
        return self.depletion_region(x)

    def get_depletion_mask(self, x=None, y=None):
        ''' Returns true for all points outside of the depletion zone
        '''
        if x is None or y is None:
            x = np.array(self.pot_data.mesh.x)
            y = np.array(self.pot_data.mesh.y)

        mask = np.zeros_like(x, dtype=np.bool)
        mask[y > self.get_depletion(x)] = True

        return mask

    def get_potential(self, x, y):
        return self.potential_grid_inter(x, y)

    def get_potential_smooth(self, x, y):
        if self.pot_smooth is None:
            self._smooth_potential()
        return self.pot_smooth(x, y, grid=False)

    def get_field(self, x, y):
        ''' Returns the field in V/um at different positions.

        Parameters
        ----------
        x, y : array_like
            Particle x, y positions
        '''

        if self.field_x is None or self.field_y is None:
            self._derive_field()
        return np.array([self.field_x(x, y, grid=False),
                         self.field_y(x, y, grid=False)])

    def _smooth_potential(self, smoothing=None):
        ''' This function takes the potential grid interpolation
            and smooths the data points.

            Smoothing is really buggy in scipy, the only
            working way is to smooth on a grid. Thus mesh points
            of the potential solution cannot be used directly.
        '''
        _LOGGER.debug('Calculate smoothed potential description')

        if not smoothing:
            smoothing = self.smoothing

        # Scale potential to make interpolation independent of bias
        v_min = np.nanmin(self.potential_grid)
        v_max = np.nanmax(self.potential_grid)

        # Scale potential to be within 0 .. 1
        potential_scaled = (self.potential_grid - v_min) / (v_max - v_min)

        def interpolator(x, y, **kwarg):
            func = RectBivariateSpline(self._x, self._y, potential_scaled.T,
                                       s=smoothing, kx=3, ky=3)
#             func = interp2d(self._x, self._y, potential_scaled, kind='cubic')
#             return func(x, y) * (v_max - v_min) + v_min
            return func(x, y, **kwarg) * (v_max - v_min) + v_min

        # Smooth on the interpolated grid
        self.pot_smooth = interpolator

    def _derive_field(self):
        ''' Takes the potential to calculate the field in x, y
        via E_x, E_y = - grad(Potential)
        with spline interpolation and smoothing.
        '''
        _LOGGER.debug('Calculate field from potential')

        if not self.pot_smooth:
            self._smooth_potential()

        E_x, E_y = np.gradient(-self.pot_smooth(self._x,
                                                self._y,
                                                grid=True),
                               np.diff(self._x)[0], np.diff(self._y)[0])

        # Create spline interpolators for E_x,E_y
        self.field_x = RectBivariateSpline(
            self._x, self._y, E_x, s=0, kx=3, ky=3)
        self.field_y = RectBivariateSpline(
            self._x, self._y, E_y, s=0, kx=3, ky=3)

    def _interpolate_nan(self, a):
        ''' Fills nans with closest non nan value.
        Might not work well for multi dimensional arrays. :TODO:
        '''
        mask = np.isnan(a)
        a[mask] = np.interp(
            np.flatnonzero(mask), np.flatnonzero(~mask), a[~mask])
        return a


def calculate_planar_sensor_w_potential(mesh, width, pitch,
                                        n_pixel, thickness):
    ''' Calculates the weighting field of a planar sensor.
    '''
    _LOGGER.info('Calculating weighting potential')
    # Mesh validity check
    mesh_width = mesh.getFaceCenters()[0, :].max(
    ) - mesh.getFaceCenters()[0, :].min()

    if mesh_width != width * n_pixel:
        raise ValueError(
            'The provided mesh width does not correspond to the sensor width')

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
        bcs.append(fipy.FixedValue(value=1.0 if pixel_position == 0. else 0.,
                                   faces=electrodes &
                                   (X > pixel_position - pitch / 2.) &
                                   (X < pixel_position + pitch / 2.)))

    solver.solve(
        potential, equation=potential.equation, boundaryConditions=bcs)
    return potential


def calculate_planar_sensor_potential(mesh, width, pitch, n_pixel, thickness,
                                      n_eff, V_bias, V_readout, V_bi=0):
    ''' Calculates the potential of a planar sensor.

        Parameters
        ----------
        mesh : fipy.Gmsh2D
               Mesh where to solve the poisson equation

        width : number
            Width of one pixel in :math:`\mathrm{\mu m}`

        pitch : number
            The width of the readout electrode in :math:`\mathrm{\mu m}`

        n_pixel : int
            Number of pixels

        thickness : number
            Thickness of the sensor in :math:`\mathrm{\mu m}`

        n_eff : number
            Effective doping concentration in :math:`\mathrm{\frac{1}{cm^3}}`

        V_bias : number
            Bias voltage in Volt

        V_readout : number
            Readout voltage in Volt

        V_bi : number
            Build in voltage. Can be calculated by
            scarce.silicon.get_diffusion_potential(...).

        Notes
        -----
        So far the depletion zone in the case of a underdepleted sensor is only
        calculated as a constant y boundary. This is wrong for pixels with low
        fill factor.
    '''

    _LOGGER.info('Calculating potential')

    # Mesh validity check
    min_x = float(mesh.getFaceCenters()[0, :].min())
    max_x = float(mesh.getFaceCenters()[0, :].max())
    min_y = float(mesh.getFaceCenters()[1, :].min())
    max_y = float(mesh.getFaceCenters()[1, :].max())

    mesh_width = max_x - min_x

    if mesh_width != width * n_pixel:
        raise ValueError(
            'The provided mesh width does not correspond to the sensor width')

    if min_y != 0:
        raise ValueError('The provided mesh does not start at 0.')

    if max_y != thickness:
        raise ValueError('The provided mesh does not end at sensor thickness.')

    # Simply add build in potential to bias potential, although that might be
    # not correct. The analytic formular does it like this.
    V_bias += V_bi

    def calculate_potential(depletion_mask=None, y_dep_new=None):
        potential = fipy.CellVariable(mesh=mesh, name='potential', value=0.)
        electrons = fipy.CellVariable(mesh=mesh, name='e-')
        electrons.valence = -1
        charge = electrons * electrons.valence
        charge.name = "charge"

        # Uniform charge distribution by setting a uniform concentration of
        # electrons = 1
        electrons.setValue(rho_scale)

        # A depletion zone within the bulk requires an internal boundary
        # condition. Internal boundary conditions seem to challenge fipy, see:
        # http://www.ctcms.nist.gov/fipy/documentation/USAGE.html#applying-internal-boundary-conditions

        large_value = 1e+10  # Hack for optimizer

        if depletion_mask is not None:
            # FIXME: Generic depletion_mask not working
            # Is overwritten here with a simple 1D depletion mask
            depletion_mask = np.logical_and(
                potential.mesh.y > y_dep_new[0],
                potential.mesh.y > y_dep_new[0])
            potential.equation = (
                fipy.DiffusionTerm(coeff=epsilon_scaled) + charge ==
                fipy.ImplicitSourceTerm(depletion_mask * large_value) -
                depletion_mask * large_value * V_bias)
        else:
            potential.equation = (
                fipy.DiffusionTerm(coeff=epsilon_scaled) + charge == 0.)

        # Calculate boundaries
        backplane = mesh.getFacesTop()
        readout_plane = mesh.getFacesBottom()

        electrodes = readout_plane
        bcs = [fipy.FixedValue(value=V_bias, faces=backplane)]
        X, _ = mesh.getFaceCenters()
        for pixel in range(n_pixel):
            pixel_position = width * (pixel + 1. / 2.) - width * n_pixel / 2.
            bcs.append(
                fipy.FixedValue(value=V_readout if pixel_position == 0.
                                else 0.,
                                faces=electrodes &
                                (X > pixel_position - pitch / 2.) &
                                (X < pixel_position + pitch / 2.)))

        solver.solve(
            potential, equation=potential.equation, boundaryConditions=bcs)
        return potential

    def get_potential(max_iter):
        r''' Solves the poisson equation for different depletion depths.

        Iterates until the minimum potential equals the bias potential.
        At this point the depletion width is correctly calculated.
        '''

        # Depletion zone info
        nx_y = 1000
        x_dep = np.linspace(min_x, max_x, nx_y)
        # Start with full depletion assumption
        y_dep = np.ones(shape=(nx_y,)) * thickness
        y_dep_new = y_dep

        description = None

        for i in range(max_iter):
            # First solution with full depletion assumption
            if i == 0:
                depletion_mask = None
            else:
                depletion_mask = description.get_depletion_mask()

            potential = calculate_potential(
                depletion_mask=depletion_mask, y_dep_new=y_dep_new)
            potential.depletion = [x_dep, y_dep]

            description = Description(potential,
                                      min_x=min_x,
                                      max_x=max_x,
                                      min_y=min_y,
                                      max_y=max_y,
                                      nx=nx_y,
                                      ny=nx_y)

#             import matplotlib.pyplot as plt
#             y = np.linspace(0, thickness, 100)
#             plt.plot(y, description.get_potential(np.zeros_like(y), y), '-', label='Pot, Numerical solution')
#             plt.plot(y, description.get_potential_smooth(0., y)[:, 0], '--', label='Pot, Numerical solution smooth')
#             plt.plot([description.get_potential_minimum_pos_y()[50], description.get_potential_minimum_pos_y()[50]], plt.ylim())
#             plt.legend(loc=0)
#             plt.show()

            potential_min = description.get_potential_minimum()
            y_dep_new = description.get_potential_minimum_pos_y()

            if ((i == 0 and np.allclose(potential_min, V_bias, rtol=1.e-2)) or
                    (i > 0 and np.allclose(potential_min, V_bias)) or
                    y_dep_new[0] >= y_dep[0]):
                potential.depletion = [description._x, y_dep]
                return potential

            # Set the new depletion depth, can only be smaller
            # Otherwise it is numerical instability
            y_dep[y_dep_new < y_dep] = y_dep_new[y_dep_new < y_dep]

        raise RuntimeError(
            'Depletion region in underdepleted sensor could not be determined')

    # The field scales with rho / epsilon, thus scale to proper value to
    # counteract numerical instabilities
    rho = constants.elementary_charge * n_eff * \
        (1e-4) ** 3  # Charge density in C / um3
    epsilon = C.epsilon_s * 1e-6  # Permitticity of silicon in F/um
    epsilon_scaled = 1.
    rho_scale = rho / epsilon

    return get_potential(max_iter=10)


def calculate_3D_sensor_potential(mesh, width_x, width_y, n_pixel_x, n_pixel_y,
                                  radius, nD, n_eff,
                                  V_bias, V_readout, V_bi=0):
    ''' Calculates the potential of a planar sensor.

        Parameters
        ----------
        mesh : fipy.Gmsh2D
               Mesh where to solve the poisson equation
        width_x : number
                  Width in x of one pixel in :math:`\mathrm{\mu m}`
        width_y : number
                  Width in y of one pixel in :math:`\mathrm{\mu m}`
        n_pixel_x : int
            Number of pixels in x
        n_pixel_y : int
            Number of pixels in y
        radius : number
                  Radius of the columns in :math:`\mathrm{\mu m}`
        nD : number
                Number of readout columns per pixel
        n_eff : number
            Effective doping concentration in :math:`\mathrm{\frac{1}{cm^3}}`
        V_bias : number
            Bias voltage in Volt
        V_readout : number
            Readout voltage in Volt
        V_bi : number
            Build in voltage. Can be calculated by
            scarce.silicon.get_diffusion_potential()

        Notes
        -----
        So far the depletion zone cannot be calculated and a fully
        depleted sensor is assumed.
    '''

    _LOGGER.info('Calculating potential')

    # Mesh validity check
    mesh_width = mesh.getFaceCenters()[0, :].max(
    ) - mesh.getFaceCenters()[0, :].min()
    mesh_height = mesh.getFaceCenters()[1, :].max(
    ) - mesh.getFaceCenters()[1, :].min()

    desc = geometry.SensorDescription3D(
        width_x, width_y, n_pixel_x, n_pixel_y, radius, nD)
    min_x, max_x, min_y, max_y = desc.get_array_corners()

    if mesh_width != max_x - min_x:
        raise ValueError(
            'Provided mesh width does not correspond to the sensor width')
    if mesh_height != max_y - min_y:
        raise ValueError(
            'Provided mesh height does not correspond to the sensor height')
    if (mesh.getFaceCenters()[0, :].min() != min_x or
            mesh.getFaceCenters()[0, :].max() != max_x):
        raise ValueError('The provided mesh has a wrong x position')
    if (mesh.getFaceCenters()[1, :].min() != min_y or
            mesh.getFaceCenters()[1, :].max() != max_y):
        raise ValueError('The provided mesh has a wrong y position')

    # The field scales with rho / epsilon, thus scale to proper value to
    # counteract numerical instabilities
    rho = constants.elementary_charge * n_eff * \
        (1e-4) ** 3  # Charge density in C / um3
    epsilon = C.epsilon_s * 1e-6  # Permitticity of silicon in F/um
    epsilon_scaled = 1.
    rho_scale = rho / epsilon

    # Define cell variables
    potential = fipy.CellVariable(mesh=mesh, name='potential', value=0.)
    electrons = fipy.CellVariable(mesh=mesh, name='e-')
    electrons.valence = -1
    charge = electrons * electrons.valence
    charge.name = "charge"

    # Uniform charge distribution by setting a uniform concentration of
    # electrons = 1
    electrons.setValue(rho_scale)

    # Add build in potential to bias potential, although that might not be
    # correct. The analytic formular does it like this
    V_bias += V_bi

    def get_potential(max_iter=10):
        ''' Calculates the potential with boundary conditions.

        Can have a larger bias column radius to simulate a not fully depleted
        sensor
        '''

        r_bias = radius  # Start with full depletion assumption

        for i in range(max_iter):
            # Set boundary condition
            bcs = []
            allfaces = mesh.getExteriorFaces()
            X, Y = mesh.getFaceCenters()

            # Set boundary conditions
            # Set readout pillars potentials
            for pos_x, pos_y in desc.get_ro_col_offsets():
                ring = allfaces & (
                    (X - pos_x) ** 2 + (Y - pos_y) ** 2 < (radius * 2) ** 2)
                bcs.append(fipy.FixedValue(value=V_readout, faces=ring))

            depletion_mask = None

            # Full bias pillars potentials = V_bias
            for pos_x, pos_y in desc.get_center_bias_col_offsets():
                ring = allfaces & (
                    (X - pos_x) ** 2 + (Y - pos_y) ** 2 < (radius) ** 2)
                bcs.append(fipy.FixedValue(value=V_bias, faces=ring))
                if not np.any(depletion_mask):
                    depletion_mask = (potential.mesh.x - pos_x)**2 + \
                        (potential.mesh.y - pos_y)**2 < r_bias ** 2
                else:
                    depletion_mask |= (potential.mesh.x - pos_x)**2 + \
                        (potential.mesh.y - pos_y)**2 < (r_bias) ** 2

            # Side bias pillars potentials = V_bias
            for pos_x, pos_y in desc.get_side_bias_col_offsets():
                ring = allfaces & (
                    (X - pos_x) ** 2 + (Y - pos_y) ** 2 < (radius) ** 2)
                bcs.append(fipy.FixedValue(value=V_bias, faces=ring))
                if not np.any(depletion_mask):
                    depletion_mask = (potential.mesh.x - pos_x)**2 + \
                        (potential.mesh.y - pos_y)**2 < r_bias ** 2
                else:
                    depletion_mask |= (potential.mesh.x - pos_x)**2 + \
                        (potential.mesh.y - pos_y)**2 < (r_bias) ** 2

            # Edge bias pillars potentials = V_bias
            for pos_x, pos_y in desc.get_edge_bias_col_offsets():
                ring = allfaces & (
                    (X - pos_x) ** 2 + (Y - pos_y) ** 2 < (radius) ** 2)
                bcs.append(fipy.FixedValue(value=V_bias, faces=ring))
                if not np.any(depletion_mask):
                    depletion_mask = (potential.mesh.x - pos_x)**2 + \
                        (potential.mesh.y - pos_y)**2 < r_bias ** 2
                else:
                    depletion_mask |= (potential.mesh.x - pos_x)**2 + \
                        (potential.mesh.y - pos_y)**2 < (r_bias) ** 2

            # A depletion zone within the bulk requires an internal boundary
            # condition. Internal boundary conditions seem to challenge fipy
            # http://www.ctcms.nist.gov/fipy/documentation/USAGE.html#applying-internal-boundary-conditions

            large_value = 1e+10  # Hack for optimizer

            potential.equation = (
                fipy.DiffusionTerm(coeff=epsilon_scaled) + charge ==
                fipy.ImplicitSourceTerm(depletion_mask * large_value) -
                depletion_mask * large_value * V_bias)

            solver.solve(potential, equation=potential.equation,
                         boundaryConditions=bcs)

            # Check if fully depleted
            if not np.isclose(potential.arithmeticFaceValue().min(),
                              V_bias, rtol=0.05, atol=0.01):
                if i == 0:
                    logging.warning('Sensor is not fully depleted. '
                                    'Try to find depletion region. ')
            else:
                return potential

            # Get line between readout and bias column to check for full
            # depletion
            for x, y in desc.get_ro_col_offsets():
                if desc.position_in_center_pixel(x, y):
                    x_ro, y_ro = x, y
                    break
            for x, y in list(desc.get_center_bias_col_offsets()) + desc.get_edge_bias_col_offsets():
                if desc.position_in_center_pixel(x, y):
                    x_bias, y_bias = x, y
                    break

            pot_descr = Description(potential,
                                    min_x=min_x,
                                    max_x=max_x,
                                    min_y=min_y,
                                    max_y=max_y,
                                    nx=width_x * n_pixel_x,
                                    ny=width_y * n_pixel_y)

            N = 1000
            x = np.linspace(x_ro, x_bias, N)
            y = np.linspace(y_ro, y_bias, N)
            # Deselect position that is within the columns
            sel = ~desc.position_in_column(x, y)
            x, y = x[sel], y[sel]
            position = np.sqrt(x ** 2 + y ** 2)  # [um]
            phi = pot_descr.get_potential(x, y)
            x_r = position.max()
            x_min = position[np.atleast_1d(np.argmin(phi))[0]]
#             import matplotlib.pyplot as plt
#             plt.plot(position, phi, color='blue', linewidth=2,
#                      label='Potential')
#             plt.plot([x_r, x_r], plt.ylim())
#             plt.plot([x_min, x_min], plt.ylim())
#             plt.show()

            # Increase bias radius boundary to simulate not depleted
            # region sourrounding bias column
            r_bias += x_r - x_min

            logging.info('Depletion region error: %d um', x_r - x_min)

        raise RuntimeError('Unable to find the depletion region')

    return get_potential(10)


def calculate_3D_sensor_w_potential(mesh, width_x, width_y,
                                    n_pixel_x, n_pixel_y, radius, nD=2):
    _LOGGER.info('Calculating weighting potential')

    # Mesh validity check
    mesh_width = mesh.getFaceCenters()[0, :].max(
    ) - mesh.getFaceCenters()[0, :].min()
    mesh_height = mesh.getFaceCenters()[1, :].max(
    ) - mesh.getFaceCenters()[1, :].min()

    desc = geometry.SensorDescription3D(
        width_x, width_y, n_pixel_x, n_pixel_y, radius, nD)
    min_x, max_x, min_y, max_y = desc.get_array_corners()

    if mesh_width != max_x - min_x:
        raise ValueError(
            'Provided mesh width does not correspond to the sensor width')
    if mesh_height != max_y - min_y:
        raise ValueError(
            'Provided mesh height does not correspond to the sensor height')
    if (mesh.getFaceCenters()[0, :].min() != min_x or
            mesh.getFaceCenters()[0, :].max() != max_x):
        raise ValueError('The provided mesh has a wrong x position')
    if (mesh.getFaceCenters()[1, :].min() != min_y or
            mesh.getFaceCenters()[1, :].max() != max_y):
        raise ValueError('The provided mesh has a wrong y position')

    potential = fipy.CellVariable(mesh=mesh, name='potential', value=0.)
    permittivity = 1.
    potential.equation = (fipy.DiffusionTerm(coeff=permittivity) == 0.)

    bcs = []
    allfaces = mesh.getExteriorFaces()
    X, Y = mesh.getFaceCenters()

    # Set boundary conditions
    # Set readout pillars potentials
    for pos_x, pos_y in desc.get_ro_col_offsets():
        ring = allfaces & ((X - pos_x) ** 2 + (Y - pos_y) ** 2 < (radius) ** 2)
        # Center pixel, phi = 1
        if desc.position_in_center_pixel(pos_x, pos_y):
            bcs.append(fipy.FixedValue(value=1., faces=ring))
        else:  # Other pixel, phi = 0
            bcs.append(fipy.FixedValue(value=0., faces=ring))

    # Full bias pillars potentials = 0
    for pos_x, pos_y in desc.get_center_bias_col_offsets():
        ring = allfaces & ((X - pos_x) ** 2 + (Y - pos_y) ** 2 < (radius) ** 2)
        bcs.append(fipy.FixedValue(value=0., faces=ring))

    # Side bias pillars potentials = 0
    for pos_x, pos_y in desc.get_side_bias_col_offsets():
        ring = allfaces & ((X - pos_x) ** 2 + (Y - pos_y) ** 2 < (radius) ** 2)
        bcs.append(fipy.FixedValue(value=0., faces=ring))

    # Edge bias pillars potentials = 0
    for pos_x, pos_y in desc.get_edge_bias_col_offsets():
        ring = allfaces & ((X - pos_x) ** 2 + (Y - pos_y) ** 2 < (radius) ** 2)
        bcs.append(fipy.FixedValue(value=0., faces=ring))

    solver.solve(
        potential, equation=potential.equation, boundaryConditions=bcs)
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
        return -1. / np.pi * (np.arctan(np.tan(ybar / 2) *
                                        np.tanh((xbar + wbar / 2.) / 2.)) -
                              np.arctan(np.tan(ybar / 2) *
                                        np.tanh((xbar - wbar / 2.) / 2.)))
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
    r''' Calculates the analytical weighting field.

        Parameters
        ----------
        x : number
            X position in sensor in um
        y : number
            Y position in sensor in um
        D : number
            The sensor thickness in um
        S : number
            3D sensor only: electrode radius in um
        is_planar : {True, False}
            To select 3D or planar geometry

        Returns
        -------
        array_like
            Weighting field in x/y in V/um

        Notes
        -----
        From Nuclear Instruments and Methods in Physics Research A 535 (2004)
        554-557, with correction

        .. math:: \bar{w} = \frac{\pi*w}{2D} -> \bar{w} = \frac{\pi*w}{D}

        The field is calculated from the drivation of the potential.
    '''

    if is_planar:
        xbar = np.pi * x / D
        ybar = np.pi * (y - D) / D
        wbar = np.pi * S / D

        # Not easy to find a more simple form
        denom = (np.cosh(1. / 2. * (wbar - 2. * xbar)) + np.cos(ybar)) * \
            (np.cosh(1. / 2. * (wbar + 2. * xbar)) + np.cos(ybar)) * D

        with np.errstate(divide='ignore', invalid='ignore'):
            E_x = - np.sin(ybar) * np.sinh(wbar / 2.) * np.sinh(xbar) / denom

            E_y = np.sinh(wbar / 2.) * \
                (np.cosh(wbar / 2.) + np.cos(ybar) * np.cosh(xbar)) / denom

        return E_x, E_y
    else:
        # 3D sensor:
        # From analytical derivation of the get_weighting_potential function
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


def get_potential_planar_analytic_1D(x, V_bias, V_readout, n_eff, D):
    r""" Calculates the potential in the depletion zone of a planar sensor.

        Parameters
        ----------
        x : array_like
            Position in the sensor between 0 and `D` in :math:`\mathrm{\mu m}`

        V_bias : number
            Bias voltage in volt.

        V_readout : number
            Readout voltage in volt.

        n_eff : number
            Effective doping concetration in :math:`\mathrm{cm^{-3}}`

        D : number
            Thickness of the sensor in :math:`\mathrm{\mu m}`

        Notes
        -----
        The formular can be derived from the 1D Poisson equation :eq:`poisson`,
        wich has the following general solution for
        :math:`x <= x_{\mathrm{dep}}` with the full depletion assumption:

        .. math:: \Phi_p = \frac{\rho}{2\epsilon} x^2 + \mathrm{const_{p,1}} x
                            + \mathrm{const_{p,2}}

        For the undepleted region :math:`x > x_{\mathrm{dep}}` there is no
        spacecharge (:math:`\rho = 0`). Thus the generel solution of the 1D
        Laplace equation :eq:`laplace` can be used here:

        .. math:: \Phi_l = \mathrm{const_{l,1}} x + \mathrm{const_{l,2}}

        For an underdepleted sensor (:math:`x_{\mathrm{dep}} <= D`) these
        boundary conditions have to be satisfied:

          1. .. math:: \Phi_p(0) = V_{\mathrm{readout}}
          2. .. math:: \Phi_l(D) = V_{\mathrm{bias}}
          3. .. math:: \Phi_p(x_{\mathrm{dep}}) = \Phi_l(x_{\mathrm{dep}})
          4. .. math:: \frac{\partial}{\partial x} \Phi_p(x_{\mathrm{dep}}) = 0
          5. .. math:: \frac{\partial}{\partial x} \Phi_l(x_{\mathrm{dep}}) = 0

        The following simultaneous equations follow:

          1. .. math:: \Phi_p = \frac{\rho}{2\epsilon} x^2
                              + \mathrm{const_{p,1}} x + V_{\mathrm{readout}}
          2. .. math:: \Phi_l = (x - D)\cdot \mathrm{const_{l,1}}
                                  + V_{\mathrm{bias}}
          3. .. math:: \frac{\rho}{2\epsilon} x_{\mathrm{dep}}^2
                        + \mathrm{const_{p,1}} x_{\mathrm{dep}}
                        + V_{\mathrm{readout}} = (x_{\mathrm{dep}} - D)
                        \cdot \mathrm{const_{l,1}} + V_{\mathrm{bias}}
          4. .. math:: \frac{\rho}{\epsilon} x_{\mathrm{dep}}
                          + \mathrm{const_{p,1}} = 0
          5. .. math:: \mathrm{const_{l,1}} = 0

        With the solution:

          .. math::

            \Phi(x) =
            \left\{
            \begin{array}{ll}
                  \frac{\rho}{2\epsilon} x^2 - \frac{\rho}{\epsilon} x_{\mathrm{dep}} x
                  + V_{\mathrm{readout}} & x\leq x_{\mathrm{dep}} \\
                  V_{\mathrm{bias}} & x > x_{\mathrm{dep}}
            \end{array}
            \right.

        with :math:`x_{\mathrm{dep}} = \sqrt{\frac{2\epsilon}{\rho}(V_{\mathrm{readout}} - V_{\mathrm{bias})}}`

        If the sensor is fully depleted (:math:`x_{\mathrm{dep}} > D`) only :math:`\Phi_p` has to be solved with the following boundary conditions:

        1. .. math:: \Phi_p(0) = V_{\mathrm{readout}}
        2. .. math:: \Phi_p(D) = V_{\mathrm{bias}}

        The following simultaneous equations follow:

        1. .. math:: \Phi_p = \frac{\rho}{2\epsilon} x^2 + \mathrm{const_{p,1}} x + V_{\mathrm{readout}}
        2. .. math:: \Phi_p(D) = \frac{\rho}{2\epsilon} D^2 + \mathrm{const_{p,1}} D + V_{\mathrm{readout}} = V_{\mathrm{bias}}

        With the solution:

          .. math::

             \Phi(x) = \frac{\rho}{2\epsilon} x^2 + \left(\frac{V_{\mathrm{bias}} - V_{\mathrm{readout}}}{D} - \frac{\rho}{2\epsilon} D\right) x + V_{\mathrm{readout}}

        For the generell solution follows:

        .. math::

            \Phi(x) =
            \left\{
            \begin{array}{ll}
                  \frac{\rho}{2\epsilon} x^2 - \frac{\rho}{\epsilon} x_{\mathrm{dep}} x + V_{\mathrm{readout}} & x_{\mathrm{dep}} \leq D, x\leq x_{\mathrm{dep}} \\
                  V_{\mathrm{bias}} & x_{\mathrm{dep}} \leq D, x > x_{\mathrm{dep}} \\
                  \frac{\rho}{2\epsilon} x^2 - \frac{\rho}{2\epsilon} D \left(\frac{x_{\mathrm{dep}}^2}{D^2} + 1\right) x + V_{\mathrm{readout}} & x_{\mathrm{dep}} > D
            \end{array}
            \right.

        with :math:`x_{\mathrm{dep}} = \sqrt{\frac{2\epsilon}{\rho}(V_{\mathrm{readout}} - V_{\mathrm{bias})}}`

    """

    # From n_eff in cm^3 to n_eff in m^3
    rho = n_eff * constants.elementary_charge * 1e-6

    x_dep = np.sqrt(2. * C.epsilon_s / rho * (V_readout - V_bias))

    # Init result array
    V = np.ones_like(x) * V_bias

    if x_dep <= D:  # Underdepleted case
        # Set spacecharge region only since not depleted area is already at
        # V_bias
        V[x <= x_dep] = rho / (2. * C.epsilon_s) * \
            x[x <= x_dep] ** 2 - rho / C.epsilon_s * x_dep * \
            x[x <= x_dep] + V_readout
    else:  # Full depleted case
        V = rho / (2. * C.epsilon_s) * x ** 2 - rho / (2. * C.epsilon_s) * \
            D * (x_dep ** 2 / D ** 2 + 1) * x + V_readout

    return V


def get_electric_field_analytic(x, y, V_bias, n_eff, D,
                                V_readout=0,
                                S=None, is_planar=True):
    r''' Calculates the electric field.

    - Planar sensor:
        Calculates the field E_y[V/um], E_x = 0 as a
        function of the position x between the electrodes [um],
        the bias voltage V_bias [V], the effective doping
        concentration n_eff [10^12 /cm^-3] and the sensor Width D [um].
        The analytical function from the detector book p. 93 is used.

        The field is derived from the potential:

        .. math:: \vec{E} = -\nabla\Phi

        From the solution of the potential follows in 1D (:math:`\vec{E}(\vec{x}) = E(x)`):

        .. math::

            E(x) =
            \left\{
            \begin{array}{ll}
                  \frac{\rho}{\epsilon} x_{\mathrm{dep}} - \frac{\rho}{\epsilon} x & x_{\mathrm{dep}} \leq D, x\leq x_{\mathrm{dep}} \\
                  0 & x_{\mathrm{dep}} \leq D, x > x_{\mathrm{dep}} \\
                  \frac{\rho}{2\epsilon} D \left(\frac{x_{\mathrm{dep}}^2}{D^2} + 1\right) - \frac{\rho}{\epsilon} x & x_{\mathrm{dep}} > D
            \end{array}
            \right.

        with :math:`x_{\mathrm{dep}} = \sqrt{\frac{2\epsilon}{\rho}(V_{\mathrm{readout}} - V_{\mathrm{bias})}}`

    - 3D sensor:
        Calculates the field E_x/E_y [V/um] in a 3d sensor as a function of the position
        x,y between the electrodes [um], the bias Voltage V_bias [V], the effective
        doping concentration n_eff [cm^-3], the electrode distance D [um] and radius R [um].
        So far the same field like the weighting field is used --> space charge is ignored.

    '''

    if is_planar:
        if S:
            raise NotImplementedError(
                'Electrode width cannot be set, \
                only full fill factor has a analytical solution!')

        # From n_eff in cm^3 to n_eff in m^3
        rho = n_eff * constants.elementary_charge * 1e-6
        x_dep = np.sqrt(2. * C.epsilon_s / rho * (V_readout - V_bias))

        # 1D case where y direction is the only dimension
        x = y

        # Init result array
        E = np.zeros_like(x)

        if x_dep <= D:  # Underdepleted case
            # Set spacecharge region only since not depleted area is already at
            # 0
            E[x <= x_dep] = rho / C.epsilon_s * x_dep - \
                rho / C.epsilon_s * x[x <= x_dep]
        else:  # Full depleted case
            E = rho / (2. * C.epsilon_s) * D * (x_dep ** 2 / D ** 2 + 1) - \
                rho / C.epsilon_s * x

        return np.zeros_like(E), E

        V_dep = silicon.get_depletion_voltage(n_eff, D)  # Depletion voltage
        a = (V_bias - V_dep) / D
        b = -2. * V_dep / (D ** 2)
        E_y = a - b * y
        return np.zeros_like(E_y), -E_y
    else:
        E_x, E_y = get_weighting_field_analytic(x, y, D, S, is_planar=False)
        E_x = E_x * V_bias
        E_y = E_y * V_bias

        return E_x, E_y
