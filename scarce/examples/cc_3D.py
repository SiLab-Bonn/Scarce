''' Example that calculates the collected charge.

    The collected charge is calculated as a function
    of the position in the sensor. The drift field
    takes irradiation into account.
'''

import matplotlib.pyplot as plt
from matplotlib import cm

from scarce import plot, silicon, analysis, sensor


def cc(bias, temperature, n_eff, t_e_trapping, t_h_trapping):
    pot_w_descr, pot_descr, geom_descr = sensor.sensor_3D(n_eff=n_eff,
                                                          V_bias=bias,
                                                          V_readout=V_readout,
                                                          temperature=temperature,
                                                          n_pixel_x=n_pixel_x,
                                                          n_pixel_y=n_pixel_y,
                                                          width_x=width_x,
                                                          width_y=width_y,
                                                          radius=radius,
                                                          nD=nD,
                                                          resolution=resolution,
                                                          smoothing=smoothing)

    return analysis.get_charge_3D(geom_descr, pot_descr, pot_w_descr,
                                  t_e_trapping=t_e_trapping, t_h_trapping=t_h_trapping,
                                  grid_x=5, grid_y=5, n_pairs=10, dt=0.001,
                                  n_steps=20000, temperature=temperature), pot_descr

if __name__ == '__main__':
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")

    n_pixel_x, n_pixel_y = 3, 3
    width_x = 250.
    width_y = 50.
    radius = 6.
    nD = 2  # Number of columns per pixel

    T = 300
    V_readout = 0.

    resolution = 80
    smoothing = 0.1

    fluence = 000.
    bias = -75
    n_eff_0 = 4.6e11
    temperature = 300

    n_eff = silicon.get_eff_acceptor_concentration(fluence, n_eff_0 / 1e12,
                                                   is_ntype=False,
                                                   is_oxygenated=True)[0] * 1e12
    t_e_trapping = silicon.get_trapping(fluence * 1e12, is_electron=True, paper=1)
    t_h_trapping = silicon.get_trapping(fluence * 1e12, is_electron=False, paper=1)

    t_e_trapping = 300000  # 4.8
    t_h_trapping = 300000  # 4.8

    (edge_x, edge_y, charge), pot_descr = cc(bias=bias,
                                             temperature=temperature,
                                             n_eff=n_eff,
                                             t_e_trapping=t_e_trapping,
                                             t_h_trapping=t_h_trapping)

    # Plot numerical potential result in 2D
    plot.plot_3D_sensor(width_x, width_y,
                        radius, nD,
                        n_pixel_x, n_pixel_y,
                        V_bias=bias, V_readout=V_readout,
                        pot_func=pot_descr.get_potential_smooth,
                        field_func=pot_descr.get_field,
                        # Comment in if you want to see the mesh
                        mesh=None,  # pot_descr.pot_data.mesh,
                        title='Potential and field of a 3D sensor, '\
                        '%dx%d pixel matrix, numerical solution' % \
                        (n_pixel_x, n_pixel_y))

    # Plot collected charge map
    plt.clf()
    plt.gca().set_aspect('equal')
    plt.gca().invert_yaxis()
    cmap = cm.get_cmap('inferno')
    cmap.set_bad('white')
    cmesh = plt.pcolormesh(edge_x, edge_y, charge,
                           cmap=cmap, vmin=0, vmax=1.05)
    plt.title('Charge collection, fluence %1.2f neq_cm2' % fluence)
    plt.grid()
    cax = plt.gcf().add_axes([plt.gca().get_position().xmax, 0.1, 0.05,
                              plt.gca().get_position().ymax - plt.gca().get_position().ymin])
    plt.colorbar(cmesh, cax=cax, orientation='vertical')
    plt.grid()
    plt.savefig('CC_%d_%d.pdf' % (fluence, bias), layout='tight')
    plt.show()
