''' Functions to create a planar_sensor or 3D sensor. '''

import logging

from scarce import (silicon, fields, geometry)

_LOGGER = logging.getLogger(__name__)


def planar_sensor(n_eff, V_bias, V_readout=0., temperature=300, n_pixel=9,
                  width=50., pitch=45., thickness=200., selection=None,
                  resolution=300., nx=None, ny=None, smoothing=0.05,
                  mesh_file='planar_mesh.msh'):
    ''' Create a planar_sensor sensor pixel array.

        Parameters
        ----------
        n_eff : number
            Effective doping concentration in :math:`\mathrm{\frac{1}{cm^3}}`
        V_bias : number
            Bias voltage in Volt
        V_readout : number
            Readout voltage in Volt
        temperature : float
            Temperature in Kelvin
        n_pixel : int
            Number of pixels
        width : number
            Width of one pixel in :math:`\mathrm{\mu m}`
        pitch : number
            Pitch (redout implant width) of one pixel in :math:`\mathrm{\mu m}`
        thickness : number
            Thickness of the sensor in :math:`\mathrm{\mu m}`
        selection : string
            Selects if the weighting potential / potentials or both are
            calculated.
            If not set: calculate weighting potential and drift potential
            If drift: calculate drift potential only
            If weighting: calculate weighting potential only
        resolution : number
            Mesh resolution. Should lead to > 200000 mesh points for exact
            results.
        nx : number
            Interpolation points in x for the potentials and fields
        ny : number
            Interpolation points in y for the potentials and fields
        smoothing : number
            Smoothing parameter for the potential. Higher number leads to
            more smooth looking potential, but be aware too much smoothing
            leads to wrong results!
        mesh_file : str
            File name of the created mesh file

        Returns
        -----
        Two scarce.fields.Description objects for the weighting potential and
        potential if no specified selection.
    '''

    # Create mesh of the sensor and stores the result
    # The created file can be viewed with any mesh viewer (e.g. gmsh)
    mesh = geometry.mesh_planar_sensor(
        n_pixel=n_pixel,
        width=width,
        thickness=thickness,
        resolution=resolution,
        filename=mesh_file)

    min_x = float(mesh.getFaceCenters()[0, :].min())
    max_x = float(mesh.getFaceCenters()[0, :].max())

    # Set um resolution grid
    if not nx:
        nx = width * n_pixel
    if not ny:
        ny = thickness

    if not selection or 'drift' in selection:
        V_bi = -silicon.get_diffusion_potential(n_eff, temperature)
        # Numerically solve the Laplace equation on the mesh
        potential = fields.calculate_planar_sensor_potential(
            mesh=mesh,
            width=width,
            pitch=pitch,
            n_pixel=n_pixel,
            thickness=thickness,
            n_eff=n_eff,
            V_bias=V_bias,
            V_readout=V_readout,
            V_bi=V_bi)
        pot_descr = fields.Description(potential,
                                       min_x=min_x,
                                       max_x=max_x,
                                       min_y=0,
                                       max_y=thickness,
                                       nx=nx,
                                       ny=ny,
                                       smoothing=smoothing)
        if selection and 'drift' in selection:
            return pot_descr

    if not selection or 'weighting' in selection:
        # Numerically solve the Poisson equation on the mesh
        w_potential = fields.calculate_planar_sensor_w_potential(
            mesh=mesh,
            width=width,
            pitch=pitch,
            n_pixel=n_pixel,
            thickness=thickness)
        pot_w_descr = fields.Description(w_potential,
                                         min_x=min_x,
                                         max_x=max_x,
                                         min_y=0,
                                         max_y=thickness,
                                         nx=nx,
                                         ny=ny,
                                         smoothing=smoothing)
        if selection and 'weighting' in selection:
            return pot_w_descr

    return pot_w_descr, pot_descr


def sensor_3D(n_eff, V_bias, V_readout=0., temperature=300, n_pixel_x=3,
              n_pixel_y=3, width_x=250., width_y=50., radius=6., nD=2,
              selection=None, resolution=80., nx=None, ny=None,
              smoothing=0.1, mesh_file='3D_mesh.msh'):
    ''' Create a 3D sensor pixel array.

        Parameters
        ----------
        n_eff : number
            Effective doping concentration in :math:`\mathrm{\frac{1}{cm^3}}`
        V_bias : number
            Bias voltage in Volt
        V_readout : number
            Readout voltage in Volt
        temperature : float
            Temperature in Kelvin
        n_pixel_x : int
            Number of pixels in x
        n_pixel_y : int
            Number of pixels in y
        width_x : number
            Width of one pixel in x in :math:`\mathrm{\mu m}`
        width_y : number
            Width of one pixel in y in :math:`\mathrm{\mu m}`
        radius : number
            Radius of readout and biasing columns in :math:`\mathrm{\mu m}`
        nD : int
            Number of readout columns per pixel
        selection : string
            Selects if the weighting potential / potentials or both are
            calculated.
            If not set: calculate weighting potential and drift potential
            If drift: calculate drift potential only
            If weighting: calculate weighting potential only
        resolution : number
            Mesh resolution. Should lead to > 300000 mesh points for exact
            results.
        nx : number
            Interpolation points in x for the potentials and fields
        ny : number
            Interpolation points in y for the potentials and fields
        smoothing : number
            Smoothing parameter for the potential. Higher number leads to
            more smooth looking potential, but be aware too much smoothing
            leads to wrong results!
        mesh_file : str
            File name of the created mesh file

        Returns
        -----
        Two scarce.fields.Description objects for the weighting potential and
        potential if not specified selection and a geometry desciption object.
    '''

    if not nx:
        nx = width_x * n_pixel_x * 4
    if not ny:
        ny = width_y * n_pixel_y * 4

    mesh = geometry.mesh_3D_sensor(width_x=width_x,
                                   width_y=width_y,
                                   n_pixel_x=n_pixel_x,
                                   n_pixel_y=n_pixel_y,
                                   radius=radius,
                                   nD=nD,
                                   resolution=resolution,
                                   filename=mesh_file)

    # Describe the 3D sensor array
    geom_descr = geometry.SensorDescription3D(width_x, width_y,
                                              n_pixel_x, n_pixel_y,
                                              radius, nD)
    min_x, max_x, min_y, max_y = geom_descr.get_array_corners()

    if not selection or 'drift' in selection:
        V_bi = -silicon.get_diffusion_potential(n_eff, temperature)
        potential = fields.calculate_3D_sensor_potential(mesh,
                                                         width_x,
                                                         width_y,
                                                         n_pixel_x,
                                                         n_pixel_y,
                                                         radius,
                                                         nD,
                                                         n_eff,
                                                         V_bias,
                                                         V_readout,
                                                         V_bi)
        pot_descr = fields.Description(potential,
                                       min_x=min_x,
                                       max_x=max_x,
                                       min_y=min_y,
                                       max_y=max_y,
                                       nx=nx,  # um res.
                                       ny=ny,  # um res.
                                       smoothing=smoothing)
        if selection and 'drift' in selection:
            return pot_descr, geom_descr

    if not selection or 'weighting' in selection:
        w_potential = fields.calculate_3D_sensor_w_potential(mesh,
                                                             width_x,
                                                             width_y,
                                                             n_pixel_x,
                                                             n_pixel_y,
                                                             radius,
                                                             nD=nD)
        pot_w_descr = fields.Description(w_potential,
                                         min_x=min_x,
                                         max_x=max_x,
                                         min_y=min_y,
                                         max_y=max_y,
                                         nx=nx,
                                         ny=ny,
                                         smoothing=smoothing
                                         )
        if selection and 'weighting' in selection:
            return pot_w_descr, geom_descr

    return pot_w_descr, pot_descr, geom_descr
