''' Example creating the weighting field and potential of a 3D pixel array.

    .. NOTE::
       The weighting potential and field is only correct if the pixel is
       surrounded by other pixels, thus `n_pixel_x = n_pixel_y = 3`
'''

import numpy as np

from scarce import (plot, sensor)


def sensor_3D():
    n_pixel_x, n_pixel_y = 3, 3
    width_x = 250.
    width_y = 50.
    radius = 6.
    nD = 2  # Number of columns per pixel
    n_eff = 1e12
    temperature = 300
    V_bias = -20.
    V_readout = 0.

    pot_descr, geom_descr = sensor.sensor_3D(n_eff=n_eff,
                                             V_bias=V_bias,
                                             V_readout=V_readout,
                                             temperature=temperature,
                                             n_pixel_x=n_pixel_x,
                                             n_pixel_y=n_pixel_y,
                                             width_x=width_x,
                                             width_y=width_y,
                                             radius=radius,
                                             nD=nD,
                                             selection='weighting',
                                             resolution=80,
                                             smoothing=0.5)

    # Plot weighting potential and field in 2D and 1d
    import matplotlib.pyplot as plt
    fig = plt.figure()
    plot.get_3D_sensor_plot(fig, width_x, width_y,
                            radius, nD,
                            n_pixel_x, n_pixel_y,
                            V_bias=1, V_readout=0,
                            pot_func=pot_descr.get_potential_smooth,
                            field_func=pot_descr.get_field,
                            # Comment in if you want to see the mesh
                            mesh=None,  # potential.mesh,
                            title='Weihting potential and field of 3D sensor,'\
                            ' %dx%d pixel matrix, numerical solution' % \
                            (n_pixel_x, n_pixel_y))

    # Get line between readout and bias column
    for x, y in geom_descr.get_ro_col_offsets():
        if geom_descr.position_in_center_pixel(x, y):
            x_ro, y_ro = x, y
            break
    for x, y in list(geom_descr.get_center_bias_col_offsets()) + geom_descr.get_edge_bias_col_offsets():
        if geom_descr.position_in_center_pixel(x, y):
            x_bias, y_bias = x, y
            break

    # Plot selected line between readout and bias column
    N = 1000
    x = np.linspace(x_ro, x_bias, N)
    y = np.linspace(y_ro, y_bias, N)
    # Deselect position that is within the columns
    sel = ~geom_descr.position_in_column(x, y)
    x, y = x[sel], y[sel]
    ax = fig.get_axes()[0]
    ax.plot(x, y, '-', color='black', linewidth=2)
    plt.show()

    # Plot weighting potential and field along selected line
    phi_smooth = pot_descr.get_potential_smooth(x, y)
    field = pot_descr.get_field(x, y)
    position = np.sqrt(x ** 2 + y ** 2)  # [um]
    plt.plot(position, phi_smooth, color='blue', linewidth=2,
             label='Potential')
    plt.legend(loc=1)
    plt.ylabel('Weighting potential [V]')
    plt.twinx(plt.gca())
    field_abs = np.sqrt(field[0] ** 2 + field[1] ** 2)
    plt.plot(position, field_abs, color='red', linewidth=2, label='Field')
    plt.grid()
    plt.legend(loc=4)
    plt.ylabel('Weighting field [V/cm]')
    plt.show()

if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    sensor_3D()
