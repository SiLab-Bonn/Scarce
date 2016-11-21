''' Example that creates a 3D pixel array with a given geometry (number of pixels in x/y, height, width, and readout columns per pixel)
    and calculates the weighting potential and fields.

    .. NOTE::
       The weighting potential and field is only correct if the pixel is surrounded by other pixels, thus
       `n_pixel_x = n_pixel_y = 3`
'''

from scarce import fields, plot, geometry


def sensor_3D():
    width_x = 250.
    width_y = 50.
    n_pixel_x, n_pixel_y = 3, 3
    radius = 6.
    resolution = 100.
    nD = 2  # Number of columns per pixel

    mesh = geometry.mesh_3D_sensor(width_x=width_x,
                                   width_y=width_y,
                                   n_pixel_x=n_pixel_x,
                                   n_pixel_y=n_pixel_y,
                                   radius=radius,
                                   nD=nD,
                                   resolution=resolution)

    potential = fields.calculate_3D_sensor_w_potential(mesh,
                                                       width_x,
                                                       width_y,
                                                       n_pixel_x,
                                                       n_pixel_y,
                                                       radius,
                                                       resolution,
                                                       nD=nD)

    # Describe the 3D sensor array
    sensor_description = geometry.SensorDescription3D(width_x, width_y, n_pixel_x, n_pixel_y, radius, nD)
    min_x, max_x, min_y, max_y = sensor_description.get_array_corners()

    # Describe the result to be able to obtain field/potential at any point in space
    field_description = fields.Description(potential,
                                           min_x=min_x,
                                           max_x=max_x,
                                           min_y=min_y,
                                           max_y=max_y,
                                           nx=width_x * n_pixel_x,
                                           ny=width_y * n_pixel_y,
                                           smoothing=0.1
                                           )

    # Plot numerical result
    plot.plot_3D_sensor(width_x, width_y,
                        radius, nD,
                        n_pixel_x, n_pixel_y,
                        V_bias=0, V_readout=1,
                        potential_function=field_description.get_potential_smooth,
                        field_function=field_description.get_field,
                        mesh=None,  # potential.mesh, # Comment in if you want to see the mesh
                        title='Weighting potential and field of a 3D sensor, 3x3 pixel matrix, numerical solution')

if __name__ == '__main__':
    sensor_3D()
