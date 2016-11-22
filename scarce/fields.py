r"""The field module calculates the potentials and fields of silicon pixel sensors.

The potential is determined numerically by solving these equations on a mesh:

.. math:: \nabla^2 \Phi = 0
   :label: laplace

.. math:: \nabla^2 \Phi = \frac{\rho}{\epsilon}
   :label: poisson

For the weighting potential equation :eq:`laplace` is solved with the boundary conditions:

.. math::

   \begin{eqnarray}
      \Phi    & = & 0 \\
      \Phi_r    & = & 1
   \end{eqnarray}

The pixel readout electrode(s) are at a potential 1 and all other equipotential pixel parts
(backside, bias columns, etc.) at 0.

For the electric potential the equation :eq:`poisson` is solved with the boundary conditions:

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
   For simple cases (e.g. planar sensor with 100% fill factor) also analytical solutions are provided.
   The analytical results are also used to benchmark the numerical results in the automated unit tests.

"""

import logging
import fipy
import numpy as np

from scipy.interpolate import RectBivariateSpline, griddata
from scipy import constants

from scarce import silicon
from scarce import geometry
from scarce import constant


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
        return [self.field_x(x, y).T, self.field_y(x, y).T]

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
    ''' Calculates the weighting field of a planar sensor.
    '''
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


def calculate_planar_sensor_potential(mesh, width, pitch, n_pixel, thickness, n_eff):
    ''' Calculates the weighting field of a planar sensor.
    '''
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
    electrons = fipy.CellVariable(mesh=mesh, name='e-')
    electrons.valence = -1
    charge = electrons * electrons.valence
    charge.name = "charge"

    # Uniform charge distribution by setting a uniform concentration of electrons = 1
    electrons.setValue(n_eff)

    permittivity = constant.epsilon_s / 100.  # to F/cm

    potential.equation = (fipy.DiffusionTerm(coeff=permittivity) + charge == 0.)

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


def calculate_3D_sensor_w_potential(mesh, width_x, width_y, n_pixel_x, n_pixel_y, radius, resolution, nD=2):
    logging.info('Calculating weighting potential')

    # Mesh validity check
    mesh_width = mesh.getFaceCenters()[0, :].max() - mesh.getFaceCenters()[0, :].min()
    mesh_height = mesh.getFaceCenters()[1, :].max() - mesh.getFaceCenters()[1, :].min()

    desc = geometry.SensorDescription3D(width_x, width_y, n_pixel_x, n_pixel_y, radius, nD)
    min_x, max_x, min_y, max_y = desc.get_array_corners()

    if mesh_width != max_x - min_x:
        raise ValueError('The provided mesh width does not correspond to the sensor width')
    if mesh_height != max_y - min_y:
        raise ValueError('The provided mesh height does not correspond to the sensor height')
    if mesh.getFaceCenters()[0, :].min() != min_x or mesh.getFaceCenters()[0, :].max() != max_x:
        raise ValueError('The provided mesh has a wrong x position')
    if mesh.getFaceCenters()[1, :].min() != min_y or mesh.getFaceCenters()[1, :].max() != max_y:
        raise ValueError('The provided mesh has a wrong y position')

    potential = fipy.CellVariable(mesh=mesh, name='potential', value=0.)
    permittivity = 1.
    potential.equation = (fipy.DiffusionTerm(coeff=permittivity) == 0.)

    bcs = []
    allfaces = mesh.getExteriorFaces()
    X, Y = mesh.getFaceCenters()

    desc = geometry.SensorDescription3D(width_x, width_y, n_pixel_x, n_pixel_y, radius, nD)

    # Set boundary conditions
    # Set readout pillars potentials
    for pos_x, pos_y in desc.get_ro_col_offsets():
        ring = allfaces & ((X - pos_x) ** 2 + (Y - pos_y) ** 2 < (radius) ** 2)
        if desc.position_in_center_pixel(pos_x, pos_y):  # Center pixel, phi = 1
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

        with np.errstate(divide='ignore', invalid='ignore'):
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


def get_potential_planar_analytic(x, V_bias, V_readout, n_eff, x_dep, D):
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
            
        x_dep : number
            Depletion zone width in :math:`\mathrm{\mu m}`. 
            
        D : number
            Thickness of the sensor in :math:`\mathrm{\mu m}`
        
        Notes
        -----
        The formular can be derived from the 1D Poisson equation :eq:`poisson`, wich has
        the following general solution for :math:`x <= x_{\mathrm{dep}}`:

        .. math:: \Phi_p = \frac{\rho}{2\epsilon} x^2 + \mathrm{const_{p,1}} x + \mathrm{const_{p,2}}
        
        For the undepleted region :math:`x > x_{\mathrm{dep}}` there is no spacecharge (:math:`\rho = 0`).
        Thus the generel solution of the 1D Laplace equation :eq:`laplace` can be used here:
        
        .. math:: \Phi_l = \mathrm{const_{l,1}} x + \mathrm{const_{l,2}}
        
        These boundary conditions have to be satisfied:

          1. .. math:: \Phi_p(0) = V_{\mathrm{readout}}
          2. .. math:: \Phi_l(D) = V_{\mathrm{bias}}
          3. .. math:: \Phi_p(x_{\mathrm{dep}}) = \Phi_l(x_{\mathrm{dep}})
          4. .. math:: \frac{\partial}{\partial x} \Phi_p(x_{\mathrm{dep}}) = \frac{\partial}{\partial x} \Phi_l(x_{\mathrm{dep}})
          
        The following simultaneous equations follow:
        
          1. .. math:: \Phi_p = \frac{\rho}{2\epsilon} x^2 + \mathrm{const_{p,1}} x + V_{\mathrm{readout}} = \frac{V_{\mathrm{dep}}}{D^2} x^2 + \mathrm{const_{p,1}} x + V_{\mathrm{readout}}
          2. .. math:: \Phi_l = (x - D)\mathrm{const_{l,1}} + V_{\mathrm{bias}}
          3. .. math:: \frac{V_{\mathrm{dep}}}{D^2} x_{\mathrm{dep}}^2 + \mathrm{const_{p,1}} x_{\mathrm{dep}} + V_{\mathrm{readout}} = (x_{\mathrm{dep}} - D)\cdot \mathrm{const_{l,1}} + V_{\mathrm{bias}}
          4. .. math:: 2 \frac{V_{\mathrm{dep}}}{D^2} x_{\mathrm{dep}} + \mathrm{const_{p,1}} = \mathrm{const_{l,1}}
          
        From 3. and 4. follows:
          
          1. .. math:: \mathrm{const_{p,1}} = \frac{V_{\mathrm{bias}} - V_{\mathrm{readout}}}{D} + \frac{x_{\mathrm{dep}}V_{\mathrm{dep}}}{D^2} \left(\frac{x_{\mathrm{dep}}}{D} - 2\right)
          2. .. math:: \mathrm{const_{l,1}} = \frac{V_{\mathrm{bias}} - V_{\mathrm{readout}}}{D} + \frac{x_{\mathrm{dep}}^2V_{\mathrm{dep}}}{D^3}

    """

    V_dep = silicon.get_depletion_voltage(n_eff, x_dep)  # Depletion voltage

    x = x[::-1]

    a = (V_bias - V_dep) / x_dep
    b = -2. * V_dep / (x_dep ** 2)
    return (a - b / 2 * x) * x


def get_electric_field_analytic(x, y, V_bias, n_eff, D, S=None, is_planar=True):
    """ Calculates the 2D electric field E_x, E_y [V/um]

    Planar sensor:
        Calculates the field E_y[V/um], E_x = 0 in a planar sensor as a
        function of the position x between the electrodes [um],
        the bias voltage V_bias [V], the effective doping
        concentration n_eff [10^12 /cm^-3] and the sensor Width D [um].
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
