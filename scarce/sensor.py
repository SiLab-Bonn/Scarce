''' Functions to create a planar or 3D sensor. '''

import logging

from scarce import (silicon, fields, geometry)

_LOGGER = logging.getLogger(__name__)


def planar(n_eff, V_bias, V_readout=0., temperature=300, n_pixel=9,
           width=50., pitch=45., thickness=200., selection=None,
           resolution=300., nx=202, ny=200, smoothing=0.1):
    ''' Create a planar sensor pixel array.

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
        resolution : number
            Mesh resolution. Should lead to > 200000 mesh points for exact
            results.
        selection : string
            Selects if the weighting potential / potentials or both are
            calculated.
            If not set: calculate weighting potential and drift potential
            If drift: calculate drift potential only
            If weighting: calculate weighting potential only

        Returns
        -----
        Two scarce.fields.Description objects for the weighting potential and
        potential if no specified selection.
    '''

    V_bi = -silicon.get_diffusion_potential(n_eff, temperature)

    # Create mesh of the sensor and stores the result
    # The created file can be viewed with any mesh viewer (e.g. gmsh)
    mesh = geometry.mesh_planar_sensor(
        n_pixel=n_pixel,
        width=width,
        thickness=thickness,
        resolution=resolution,
        filename='planar_mesh.msh')

    min_x = float(mesh.getFaceCenters()[0, :].min())
    max_x = float(mesh.getFaceCenters()[0, :].max())

    if not selection or 'drift' in selection:
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
