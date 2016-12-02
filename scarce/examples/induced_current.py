''' Example that creates a planar silicon sensor with a given geometry (thickness, number of pixels, pitch, width).
    Calculates the electrical potential and fields. For comparison also the analytical result of a planar sensor with
    100% fill factor (width = pitch) is created.

    .. WARNING::
       The calculation of the depletion region is simplified. If the depletion is not at a contant y position
       in the sensor (e.g. for pixels with small fill factor) it deviates from the correct solution.
'''

import numpy as np
import matplotlib.pyplot as plt
from scarce import fields, plot, geometry, silicon, tools


def induced_current():
    # Number of pixels influences how correct the field for the
    # center pixel(s) is due to more far away infinite boundary condition
    n_pixel = 9

    # Geometry of one pixel
    width = 50.
    thickness = 300.
    pitch = 45.

    n_eff = 5e12  # n_eff [cm^-3]
    temperature = 300

    # Potentials
    V_bias = -160.
    V_readout = 0.
    V_bi = -silicon.get_diffusion_potential(n_eff, temperature)

# Create mesh of the sensor and stores the result
# The created file can be viewed with any mesh viewer (e.g. gmsh)
#     mesh = geometry.mesh_planar_sensor(
#         n_pixel=n_pixel,
#         width=width,
#         thickness=thickness,
#         resolution=200.,
#         filename='planar_mesh_example.msh')
#
# Numerically solve the laplace equation on the mesh
#     potential = fields.calculate_planar_sensor_potential(mesh=mesh,
#                                                          width=width,
#                                                          pitch=pitch,
#                                                          n_pixel=n_pixel,
#                                                          thickness=thickness,
#                                                          n_eff=n_eff,
#                                                          V_bias=V_bias,
#                                                          V_readout=V_readout,
#                                                          V_bi=V_bi)
#
#     w_potential = fields.calculate_planar_sensor_w_potential(mesh=mesh,
#                                                              width=width,
#                                                              pitch=pitch,
#                                                              n_pixel=n_pixel,
#                                                              thickness=thickness)
#
#     min_x = float(mesh.getFaceCenters()[0, :].min())
#     max_x = float(mesh.getFaceCenters()[0, :].max())
#     pot_descr = fields.Description(potential,
#                                    min_x=min_x,
#                                    max_x=max_x,
#                                    min_y=0,
#                                    max_y=thickness)
#     pot_w_descr = fields.Description(w_potential,
#                                      min_x=min_x,
#                                      max_x=max_x,
#                                      min_y=0,
#                                      max_y=thickness)
#
#     pot_descr.get_field(0, 0)
#     pot_w_descr.get_field(0, 0)
#
#     tools.save(pot_descr, 'planar_sensor_pot.sc')
#     tools.save(pot_w_descr, 'planar_sensor_pot_w.sc')

    pot_descr = tools.load('planar_sensor_pot.sc')
    pot_w_descr = tools.load('planar_sensor_pot_w.sc')

    # Get analytical result
    def potential_analytic(x, y):
        return fields.get_weighting_potential_analytic(x, y, D=thickness, S=width, is_planar=True)

#     # Plot analytical / numerical result with depletion region in 1D
#     y = np.linspace(0, thickness, 100)
#     x = np.zeros_like(y)
#     plt.plot(y, pot_descr.get_potential(x, y), label='Potential, numerical', linewidth=2)
#     pot_masked = np.ma.masked_array(pot_descr.get_potential(x, y), mask=pot_descr.get_depletion_mask(x, y), linewidth=2)
#     plt.plot(y, pot_masked, label='Potential, numerical, depleted', linewidth=2)
#     plt.plot([pot_descr.get_depletion(x[50]), pot_descr.get_depletion(x[50])], plt.ylim(), label='Depletion, numerical ', linewidth=2)
#     plt.plot(y, fields.get_potential_planar_analytic_1D(y, V_bias=V_bias + V_bi, V_readout=V_readout, n_eff=n_eff, D=thickness), '--', label='Potential, analytical', linewidth=2)
#     plt.plot([silicon.get_depletion_depth(np.abs(V_bias), n_eff / 1e12, temperature), silicon.get_depletion_depth(np.abs(V_bias), n_eff / 1e12, temperature)], plt.ylim(), '--', label='Depletion, analytical', linewidth=2)
#     plt.legend(loc=0)
#     plt.title('Potential in a not fully depleted planar sensor')
#     plt.xlabel('Position [um]')
#     plt.xlabel('Potential [V]')
#     plt.grid()
#     plt.show()

#     # Plot numerical result in 2D
#     plot.plot_planar_sensor(width=width,
#                             pitch=pitch,
#                             thickness=thickness,
#                             n_pixel=n_pixel,
#                             V_backplane=0.,
#                             V_readout=1.,
#                             potential_function=pot_w_descr.get_potential,
#                             field_function=pot_w_descr.get_field,
#                             title='Planar sensor weighting potential')
#     plot.plot_planar_sensor(width=width,
#                             pitch=pitch,
#                             thickness=thickness,
#                             n_pixel=n_pixel,
#                             V_backplane=V_bias,  # Weighting field = 0 at backplane
#                             V_readout=V_readout,  # Weighting field = 1 at readout
#                             potential_function=pot_descr.get_potential,
#                             field_function=pot_descr.get_field,
#                             depletion_function=pot_descr.get_depletion,
#                             title='Planar sensor potential')

    fields.

if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    induced_current()
