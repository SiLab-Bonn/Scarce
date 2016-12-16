''' Example that creates a 3D pixel array with a given geometry (number of pixels in x/y, height, width, and readout columns per pixel)
    and calculates the weighting potential and fields.

    .. NOTE::
       The weighting potential and field is only correct if the pixel is surrounded by other pixels, thus
       `n_pixel_x = n_pixel_y = 3`
'''

import numpy as np

from scarce import fields, plot, geometry, tools


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
                                                       nD=nD)
 
    # Describe the 3D sensor array
    sensor_descr = geometry.SensorDescription3D(width_x, width_y, n_pixel_x, n_pixel_y, radius, nD)
    min_x, max_x, min_y, max_y = sensor_descr.get_array_corners()

 
    # Describe the result to be able to obtain field/potential at any point in space
    pot_descr = fields.Description(potential,
                                           min_x=min_x,
                                           max_x=max_x,
                                           min_y=min_y,
                                           max_y=max_y,
                                           nx=width_x * n_pixel_x * 2.,
                                           ny=width_y * n_pixel_y * 2.,
                                           smoothing=0.1
                                           )
     
#     tools.save(pot_descr, 'tmp.sc')
#     pot_descr = tools.load('tmp.sc')

#     # Plot numerical result
#     plot.plot_3D_sensor(width_x, width_y,
#                         radius, nD,
#                         n_pixel_x, n_pixel_y,
#                         V_bias=0, V_readout=1,
#                         potential_function=pot_descr.get_potential_smooth,
#                         field_function=pot_descr.get_field,
#                         mesh=None,  # potential.mesh, # Comment in if you want to see the mesh
#                         title='Weighting potential and field of a 3D sensor, %dx%d pixel matrix, numerical solution' % (n_pixel_x, n_pixel_y))

    import matplotlib.pyplot as plt
#     fig = plt.figure()
#     plot.get_3D_sensor_plot(fig, width_x, width_y,
#                             radius, nD,
#                             n_pixel_x, n_pixel_y,
#                             V_bias=1, V_readout=0,
#                             potential_function=pot_descr.get_potential_smooth,
#                             field_function=pot_descr.get_field,
#                             # Comment in if you want to see the mesh
#                             mesh=None,  # potential.mesh,
#                             title='Potential and field of a 3D sensor, '\
#                             '%dx%d pixel matrix, numerical solution' % \
#                             (n_pixel_x, n_pixel_y))
#     
    for x, y in sensor_descr.get_ro_col_offsets():
        if sensor_descr.position_in_center_pixel(x, y):
            x_ro, y_ro = x, y
            break
 
    for x, y in list(sensor_descr.get_center_bias_col_offsets()) + sensor_descr.get_edge_bias_col_offsets():
        if sensor_descr.position_in_center_pixel(x, y):
            x_bias, y_bias = x, y
            break
     
    N = 1000
   
    x = np.linspace(x_ro, x_bias, N)
    y = np.linspace(y_ro, y_bias, N)

#     ax = fig.get_axes()[0] 
#     ax.plot(x, y,'-', color='black', linewidth=2)
#     
#     plt.show()
    
    sel = ~sensor_descr.position_in_column(x, y)
    
    print pot_descr.potential_grid.shape
    
    xx, yy = np.meshgrid(pot_descr._x, pot_descr._y, sparse=False)
    
    pot_descr.potential_grid[pot_descr.potential_grid > 0.999] = 1.
#     print pot_descr._xx.shape
#     print pot_descr._yy.shape
#     raise
    phi = pot_descr.get_potential(x, y)
    phi_smooth = pot_descr.get_potential_smooth(x, y)
    field = pot_descr.get_field(x, y)
    fa = np.sqrt(field[0] ** 2 + field[1] ** 2)
    plt.plot(np.arange(N), phi, label='Pot')
    plt.plot(np.arange(N)[sel], phi[sel], label='Pot sel')
    plt.plot(np.arange(N), phi_smooth, label='Pot smooth')
    plt.plot(np.arange(N), fa, label='Field')
#     plt.legend(loc=1)
#     plt.twinx(plt.gca())
#     plt.plot(np.arange(N), x)
#     plt.plot(np.arange(N), y)
#     plt.plot(np.arange(N)[sel], phi[sel], '.', label='Pot')
    
    plt.legend(loc=4)
    plt.show()

if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    sensor_3D()
