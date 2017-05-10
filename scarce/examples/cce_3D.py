''' Example that calculates the collected charge.

    The collected charge is calculated as a function
    of the position in the sensor. The drift field
    takes irradiation into account.
'''

import matplotlib.pyplot as plt
from matplotlib import cm

from scarce import plot, silicon, analysis, sensor


def cc(bias, temperature, n_eff, t_e_trapping, t_h_trapping, geom_descr, pot_descr, pot_w_descr):
    return analysis.get_charge_3D(geom_descr, pot_descr, pot_w_descr,
                                  t_e_trapping=t_e_trapping, t_h_trapping=t_h_trapping,
                                  grid_x=5, grid_y=5, n_pairs=10, dt=0.001,
                                  n_steps=20000, temperature=temperature)


def cc_efficiency(biases, n_effs, t_e_trappings, t_h_trappings, t_e_1_trapping, t_h_1_trapping):
    n_pixel_x, n_pixel_y = 3, 3
    width_x = 250.
    width_y = 50.
    radius = 6.
    nD = 2  # Number of columns per pixel

    V_readout = 0.

    resolution = 80
    smoothing = 0.1

    pot_w_descr, pot_descr, geom_descr = sensor.sensor_3D(n_eff=n_eff_0,
                                                          V_bias=bias_0,
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

    edge_x, edge_y, charge = cc(bias=-20,
                                temperature=temperature,
                                n_eff=n_eff_0,
                                t_e_trapping=None,
                                t_h_trapping=None,
                                pot_w_descr=pot_w_descr,
                                pot_descr=pot_descr,
                                geom_descr=geom_descr)

    ccs = [charge]

    # Plot collected charge map
    plt.clf()
    plt.gca().set_aspect('equal')
    plt.gca().invert_yaxis()
    cmap = cm.get_cmap('inferno')
    cmap.set_bad('white')
    cmesh = plt.pcolormesh(edge_x, edge_y, ccs[0],
                           cmap=cmap, vmin=0, vmax=1.05)
    plt.title('Charge collection, fluence %1.2f neq_cm2' % 0)
    plt.grid()
    cax = plt.gcf().add_axes([plt.gca().get_position().xmax, 0.1, 0.05,
                              plt.gca().get_position().ymax - plt.gca().get_position().ymin])
    plt.colorbar(cmesh, cax=cax, orientation='vertical')
    plt.grid()
    plt.savefig('CC_%d_%d.pdf' % (0, bias_0), layout='tight')
    plt.show()

    with open("data.txt", "a+") as myfile:
        myfile.append('Bias\tN_eff\tN_eff\tt_e\tt_h\n')

        for i, n_eff in enumerate(n_effs):
            cc = cc_efficiency(biase=biases[i],
                               n_eff=n_eff,
                               t_e_trapping=t_e_trappings[i],
                               t_h_trapping=t_h_trappings[i])

            # Plot collected charge map
            plt.clf()
            plt.gca().set_aspect('equal')
            plt.gca().invert_yaxis()
            cmap = cm.get_cmap('inferno')
            cmap.set_bad('white')
            cmesh = plt.pcolormesh(edge_x, edge_y, ccs[0],
                                   cmap=cmap, vmin=0, vmax=1.05)
            plt.title('Charge collection, fluence %1.2f neq_cm2' % 0)
            plt.grid()
            cax = plt.gcf().add_axes([plt.gca().get_position().xmax, 0.1, 0.05,
                                      plt.gca().get_position().ymax - plt.gca().get_position().ymin])
            plt.colorbar(cmesh, cax=cax, orientation='vertical')
            plt.grid()
            plt.savefig('CC_%d_%d.pdf' % (n_eff, biases[i]), layout='tight')
            plt.show()

        for i, n_eff in enumerate(n_effs):
            myfile.append('%d\t%1.2e\t%1.2e\t%1.2e\t%1.2e\n', (biases[i],
                                                               n_eff))

if __name__ == '__main__':
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")

    n_eff_0 = 4.6e11
    bias_0 = -20
    temperature = 300

    biases = [-10]
    n_effs = silicon.get_eff_acceptor_concentration(fluence=1000., n_eff_0=n_eff_0 / 1e12,
                                                    is_ntype=False,
                                                    is_oxygenated=True)[0] * 1e12
    t_e_trappings = silicon.get_trapping(1000 * 1e12, is_electron=True,
                                         paper=1)
    t_h_trappings = silicon.get_trapping(1000 * 1e12, is_electron=False,
                                         paper=1)
    t_e_1_trapping = 0.
    t_h_1_trapping = 0.

    cc_efficiency(biases, n_effs, t_e_trappings, t_h_trappings, t_e_1_trapping, t_h_1_trapping)
