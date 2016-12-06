''' Example that creates a planar silicon sensor with a given geometry.

    Calculates the weighting potential and fields. For comparison also the
    analytical result of a planar sensor with 100% fill factor (width = pitch)
    is created.

    .. NOTE::
       With increasing distance from the center pixel the numerical result
       deviates from the analytical one. This shows that is is important to
       use several pixels (> 5) to get a proper field description in the
       center pixel.
'''

import numpy as np
import matplotlib.pyplot as plt

from scarce import fields, plot, geometry


def sensor_planar():
    # Number of pixels influences how correct the field for the
    # center pixel(s) is due to more far away infinite boundary condition
    n_pixel = 9

    # Geometry of one pixel
    width = 50.
    thickness = 200.
    pitch = 45.

    # Create mesh of the sensor and stores the result
    # The created file can be viewed with any mesh viewer (e.g. gmsh)
    mesh = geometry.mesh_planar_sensor(
        n_pixel=n_pixel,
        width=width,
        thickness=thickness,
        resolution=300.,
        filename='planar_mesh_example.msh')

    # Numerically solve the laplace equation on the mesh
    potential = fields.calculate_planar_sensor_w_potential(mesh=mesh,
                                                           width=width,
                                                           pitch=pitch,
                                                           n_pixel=n_pixel,
                                                           thickness=thickness)

    # Describe result to obtain field/potential at any point in space
    desc = fields.Description(potential,
                              min_x=-width * float(n_pixel),
                              max_x=width * float(n_pixel),
                              min_y=0,
                              max_y=thickness,
                              nx=200 * n_pixel,
                              ny=200,
                              smoothing=0.1
                              )

    # Analytical results
    def potential_analytic(x, y):
        return fields.get_weighting_potential_analytic(x, y,
                                                       D=thickness, S=width,
                                                       is_planar=True)

    def field_analytic(x, y):
        return fields.get_weighting_field_analytic(x, y, D=thickness, S=width,
                                                   is_planar=True)

    # Plot analytical / numerical result with in 1D
    y = np.linspace(0, thickness, 1000)
    x = np.zeros_like(y)
    plt.plot(y, desc.get_potential_smooth(x, y),
             label='Potential, numerical', linewidth=2)
    plt.plot(y, potential_analytic(x, y),
             '--', label='Potential, analytical', linewidth=2)
    plt.ylabel('Potential [V]')
    plt.legend(loc=1)
    ax2 = plt.gca().twinx()
    ax2.plot(y, desc.get_field(x, y)[1],
             '--', label='Field, numerical', linewidth=2)
    ax2.plot(y, field_analytic(x, y)[1],
             '--', label='Field, analytical', linewidth=2)
    plt.ylabel('Field [V/$\mu m$]')
    plt.legend(loc=4)
    plt.title('Weighting potential and field in a planar sensor')
    plt.xlabel('Position y [um]')
    plt.grid()
    plt.show()

    # Plot 2D results

    # Plot numerical result
    plot.plot_planar_sensor(width=width,
                            pitch=pitch,
                            thickness=thickness,
                            n_pixel=n_pixel,
                            V_backplane=0,  # Weighting field = 0 at backplane
                            V_readout=1,  # Weighting field = 1 at readout
                            potential_function=desc.get_potential_smooth,
                            field_function=desc.get_field,
                            # Comment in if you want to see the mesh
                            mesh=None,  # potential.mesh,
                            title='Planar sensor mesh')

    # Plot analytic result
    plot.plot_planar_sensor(width=width,
                            # Analytical solution exist only for pitch =width
                            # (100% fill factor)
                            pitch=width,
                            thickness=thickness,
                            n_pixel=n_pixel,
                            V_backplane=0,  # Weighting field = 0 at backplane
                            V_readout=1,  # Weighting field = 1 at readout
                            potential_function=potential_analytic,
                            title='Planar sensor weighting potential and field, \
                                    analytical solution')

if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    sensor_planar()
