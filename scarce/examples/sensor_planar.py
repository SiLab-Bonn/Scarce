''' Example that creates a planar silicon sensor with a given geometry (thickness, number of pixels, pitch, width).
    Calculates the electrical potential and fields. For comparison also the analytical result of a planar sensor with
    100% fill factor (width = pitch) is created.

    .. WARNING::
       The calculation of the depletion region is simplified. If the depletion is not at a contant y position
       in the sensor (e.g. for pixels with small fill factor) it deviates from the correct solution.
'''

import numpy as np
import matplotlib.pyplot as plt
from scarce import fields, plot, geometry, silicon


def sensor_planar():
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
                                                         V_readout=V_readout,
                                                         V_bi=V_bi)

    min_x = float(mesh.getFaceCenters()[0, :].min())
    max_x = float(mesh.getFaceCenters()[0, :].max())
    description = fields.Description(potential,
                                     min_x=min_x,
                                     max_x=max_x,
                                     min_y=0,
                                     max_y=thickness)

    # Plot analytical / numerical result with depletion region in 1D
    y = np.linspace(0, thickness, 100)
    x = np.zeros_like(y)
    plt.plot(y, description.get_potential(x, y), label='Potential, numerical', linewidth=2)
    pot_masked = np.ma.masked_array(description.get_potential(x, y), mask=description.get_depletion_mask(x, y), linewidth=2)
    plt.plot(y, pot_masked, label='Potential, numerical, depleted', linewidth=2)
    plt.plot([description.get_depletion(x[50]), description.get_depletion(x[50])], plt.ylim(), label='Depletion, numerical ', linewidth=2)
    plt.plot(y, fields.get_potential_planar_analytic_1D(y, V_bias=V_bias + V_bi, V_readout=V_readout, n_eff=n_eff, D=thickness), '--', label='Potential, analytical', linewidth=2)
    plt.plot([silicon.get_depletion_depth(np.abs(V_bias), n_eff / 1e12, temperature), silicon.get_depletion_depth(np.abs(V_bias), n_eff / 1e12, temperature)], plt.ylim(), '--', label='Depletion, analytical', linewidth=2)
    plt.legend(loc=0)
    plt.title('Potential in a not fully depleted planar sensor')
    plt.xlabel('Position [um]')
    plt.xlabel('Potential [V]')
    plt.grid()
    plt.show()

    # Plot numerical result in 2D
    plot.plot_planar_sensor(width=width,
                            pitch=pitch,
                            thickness=thickness,
                            n_pixel=n_pixel,
                            V_backplane=V_bias,  # Weighting field = 0 at backplane
                            V_readout=V_readout,  # Weighting field = 1 at readout
                            potential_function=description.get_potential,
                            field_function=description.get_field,
                            depletion_function=description.get_depletion,
                            mesh=None,  # potential.mesh,  # Comment in if you want to see the mesh
                            title='Planar sensor potential')

if __name__ == '__main__':
    sensor_planar()
