''' Example that creates a 3D pixel array with a given geometry (number of pixels in x/y, height, width, and readout columns per pixel)
    and calculates the potential and fields.

    .. WARNING::
       The calculation of a partially depleted sensor is not supported right now.
'''

from scarce import fields, plot, geometry, silicon


def sensor_3D():
    # Number of pixels influences how correct the field for the
    # center pixel(s) is due to more correct boundary condition
    n_pixel_x, n_pixel_y = 1, 1

    # Geometry of one pixel
    width_x = 250.
    width_y = 50.
    radius = 6.
    nD = 2  # Number of columns per pixel

    n_eff = 1e12  # n_eff [cm^-3]
    temperature = 300

    # Potentials
    V_bias = -20.
    V_readout = 0.
    V_bi = -silicon.get_diffusion_potential(n_eff, temperature)

    mesh = geometry.mesh_3D_sensor(width_x=width_x,
                                   width_y=width_y,
                                   n_pixel_x=n_pixel_x,
                                   n_pixel_y=n_pixel_y,
                                   radius=radius,
                                   nD=nD,
                                   resolution=50.)

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

    # Describe the 3D sensor array
    sensor_description = geometry.SensorDescription3D(width_x, width_y, n_pixel_x, n_pixel_y, radius, nD)
    min_x, max_x, min_y, max_y = sensor_description.get_array_corners()

    # Describe the result to be able to obtain field/potential at any point in space
    field_description = fields.Description(potential,
                                           min_x=min_x,
                                           max_x=max_x,
                                           min_y=min_y,
                                           max_y=max_y,
                                           nx=width_x * n_pixel_x,  # um resolution
                                           ny=width_y * n_pixel_y,  # um resolution
                                           smoothing=0.1
                                           )

    # Plot numerical result
    plot.plot_3D_sensor(width_x, width_y,
                        radius, nD,
                        n_pixel_x, n_pixel_y,
                        V_bias=V_bias + V_bi, V_readout=V_readout,
                        potential_function=field_description.get_potential_smooth,
                        field_function=field_description.get_field,
                        mesh=None,  # potential.mesh, # Comment in if you want to see the mesh
                        title='Potential and field of a 3D sensor, %dx%d pixel matrix, numerical solution' % (n_pixel_x, n_pixel_y))

if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    sensor_3D()
