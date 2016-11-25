''' Example that creates a planar silicon sensor with a given geometry (thickness, number of pixels, pitch, width).
    Calculates the electrical potential and fields. For comparison also the analytical result of a planar sensor with
    100% fill factor (width = pitch) is created.

    .. NOTE::
       With increasing distance from the center pixel the numerical result deviates from the analytical one.
       This shows that is is important to use several pixels (> 5) to get a proper field description in the center pixel.
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
    pitch = width

    n_eff = 1e11  # n_eff [cm^-3]

    V_bias = -0.5
    V_readout = 0.

    # Create mesh of the sensor and stores the result
    # The created file can be viewed with any mesh viewer (e.g. gmsh)
    mesh = geometry.mesh_planar_sensor(
        n_pixel=n_pixel,
        width=width,
        thickness=thickness,
        resolution=200.,
        filename='planar_mesh_example.msh')

    # Numerically solve the laplace equation on the mesh
    potential = fields.calculate_planar_sensor_potential(mesh=mesh,
                                                         width=width,
                                                         pitch=pitch,
                                                         n_pixel=n_pixel,
                                                         thickness=thickness,
                                                         n_eff=n_eff,
                                                         V_bias=V_bias,
                                                         V_readout=V_readout)

    # Describe the result to be able to obtain field/potential at any point in space
    description = fields.Description(potential,
                                     min_x=-width * float(n_pixel),
                                     max_x=width * float(n_pixel),
                                     min_y=0,
                                     max_y=thickness,
                                     nx=200 * n_pixel,
                                     ny=200,
                                     smoothing=0.1
                                     )

# Plot numerical result
#     plot.plot_planar_sensor(width=width,
#                             pitch=pitch,
#                             thickness=thickness,
#                             n_pixel=n_pixel,
# V_backplane=0,  # Weighting field = 0 at backplane
# V_readout=1,  # Weighting field = 1 at readout
#                             potential_function=description.get_potential_smooth,
#                             field_function=description.get_field,
# mesh=potential.mesh,  # Comment in if you want to see the mesh
#                             title='Planar sensor mesh')

    # Get analytical result
    def potential_analytic(x, y):
        return fields.get_weighting_potential_analytic(x, y, D=thickness, S=width, is_planar=True)

    # Plot analytical / numerical result
    y = np.linspace(0, thickness, 100)
    x = np.zeros_like(y)

    xx, yy = np.meshgrid(x, y, sparse=True)
    pot_numeric = description.get_potential(xx, yy)

#     print description.get_potential_smooth(xx, yy)[50, :]

    min_x = float(mesh.getFaceCenters()[0, :].min())
    max_x = float(mesh.getFaceCenters()[0, :].max())
    description = fields.Description(potential,
                                      min_x=min_x,
                                      max_x=max_x,
                                      min_y=0,
                                      max_y=thickness)
            
    print description.get_potential_minimum(axis=0)
    print description.get_potential_minimum(axis=0).shape
    
    print description.get_potential_minimum_pos_y().shape
    
    print description.get_potential_minimum_pos_y()

    plt.plot(y, fields.get_potential_planar_analytic_1D(y, V_bias=V_bias, V_readout=V_readout, n_eff=n_eff, D=thickness), label='Pot, analytical poisson solution')
#     plt.plot(y, fields.get_electric_field_analytic(x, y, V_bias=V_bias, n_eff=n_eff, D=thickness)[1], label='Ey, analytical poisson solution')
    plt.plot(y, description.get_potential(x, y), label='Pot, Numerical solution')
    plt.plot(y, description.get_potential_smooth(0., y)[:, 0], '--', label='Pot, Numerical solution smooth')
    plt.plot([description.get_potential_minimum_pos_y()[50], description.get_potential_minimum_pos_y()[50]], plt.ylim()) 
    plt.legend(loc=0)
    plt.show()

    # Plot numerical result
    plot.plot_planar_sensor(width=width,
                            pitch=pitch,
                            thickness=thickness,
                            n_pixel=n_pixel,
                            V_backplane=V_bias,  # Weighting field = 0 at backplane
                            V_readout=V_readout,  # Weighting field = 1 at readout
                            potential_function=description.get_potential,
                            field_function=description.get_field,
                            mesh=None,  # potential.mesh,  # Comment in if you want to see the mesh
                            title='Planar sensor mesh')

#     plot.plot_planar_sensor(width=width,
# pitch=width,  # Analytical solution exist only for pitch =width (100% fill factor)
#                             thickness=thickness,
#                             n_pixel=n_pixel,
# V_backplane=0,  # Weighting field = 0 at backplane
# V_readout=1,  # Weighting field = 1 at readout
#                             potential_function=potential_analytic,
#                             title='Planar sensor weighting potential and field, analytical solution')

if __name__ == '__main__':
    sensor_planar()
