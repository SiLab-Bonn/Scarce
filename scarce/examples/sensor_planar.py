''' Example that creates a planar silicon sensor with a given geometry.

    Calculates the electrical potential and fields. For comparison also the
    analytical result of a planar sensor with 100% fill factor (width = pitch)
    is created.

    .. WARNING::
       The calculation of the depletion region is simplified. If the depletion
       is not at a contant y position in the sensor (e.g. for pixels with very
       small fill factor) it deviates from the correct solution.
'''

import numpy as np
import matplotlib.pyplot as plt
from scarce import fields, plot, silicon, sensor


def sensor_planar():
    # Sensor parameters
    n_eff = 6.2e12
    n_pixel = 9
    width = 50.
    pitch = 30.
    thickness = 250.
    smoothing = 0.05
    resolution = 287
    temperature = 300.
    V_bias = -300.
    V_readout = 0.

    # Create sensor
    pot_descr = sensor.planar_sensor(n_eff=n_eff,
                                     V_bias=V_bias,
                                     V_readout=V_readout,
                                     temperature=temperature,
                                     n_pixel=n_pixel,
                                     width=width,
                                     pitch=pitch,
                                     thickness=thickness,
                                     # Calculate drift potential only
                                     # to safe time
                                     selection='drift',
                                     resolution=resolution,
                                     # Might have to be adjusted when changing
                                     # the geometry
                                     smoothing=smoothing
                                     )

    # Build in voltage needed for analytical solution
    V_bi = -silicon.get_diffusion_potential(n_eff, temperature)

    # Plot analytical / numerical result with depletion region in 1D
    y = np.linspace(0, thickness, 1000)
    x = np.zeros_like(y)

    plt.plot(y, pot_descr.get_potential(x, y),
             label='Potential, numerical', linewidth=2)
    pot_masked = np.ma.masked_array(pot_descr.get_potential(x, y),
                                    mask=pot_descr.get_depl_mask(x, y))
    plt.plot(y, pot_masked, label='Potential, numerical, depl.',
             linewidth=2)
    plt.plot([pot_descr.get_depletion(x[500]),
              pot_descr.get_depletion(x[500])],
             plt.ylim(), label='Depletion, numerical ', linewidth=2)
    plt.plot(y, fields.get_potential_planar_analytic_1D(y,
                                                        V_bias=V_bias + V_bi,
                                                        V_readout=V_readout,
                                                        n_eff=n_eff,
                                                        D=thickness),
             '--', label='Potential, analytical', linewidth=2)
    plt.plot([silicon.get_depletion_depth(np.abs(V_bias), n_eff / 1e12,
                                          temperature),
              silicon.get_depletion_depth(np.abs(V_bias), n_eff / 1e12,
                                          temperature)],
             plt.ylim(), '--', label='Depletion, analytical', linewidth=2)
    plt.ylabel('Potential [V]')
    plt.legend(loc=1)
    plt.ylabel('Potential [V]')
    ax2 = plt.gca().twinx()
    ax2.plot(y, pot_descr.get_field(x, y)[1],
             '--', label='Field, numerical', linewidth=2)
    ax2.plot(y, fields.get_electric_field_analytic(x,
                                                   y,
                                                   V_bias=V_bias,
                                                   V_readout=V_readout,
                                                   n_eff=n_eff,
                                                   D=thickness)[1],
             '--', label='Field, analytical', linewidth=2)
    plt.ylabel('Field [V/cm]')
    plt.legend(loc=4)
    plt.ylabel('Field [V/cm]')
    plt.title('Potential in a not fully depleted planar sensor')
    plt.xlabel('Position [um]')

    plt.grid()
    plt.show()

    # Plot numerical result in 2D
    plot.plot_planar_sensor(width=width,
                            pitch=pitch,
                            thickness=thickness,
                            n_pixel=n_pixel,
                            # Weighting field = 0 at backplane
                            V_backplane=V_bias,
                            # Weighting field = 1 at readout
                            V_readout=V_readout,
                            pot_func=pot_descr.get_potential,
                            field_func=pot_descr.get_field,
                            depl_func=pot_descr.get_depletion,
                            # Comment in if you want to see the mesh
                            mesh=None,  # potential.mesh,
                            title='Planar sensor potential')

if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    sensor_planar()
