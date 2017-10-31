''' Example that calculates the collected charge.

    The collected charge is calculated as a function
    of the position in the sensor. The drift field
    takes irradiation into account.
'''

import matplotlib.pyplot as plt
from matplotlib import cm

from scarce import tools
from scarce.examples import cc_3D


if __name__ == '__main__':
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")

    n_eff_0 = 4.6e11
    n_pixel_x, n_pixel_y = 3, 3
    width_x = 250.
    width_y = 50.
    radius = 6.
    nD = 2  # Number of columns per pixels

    resolution = 80
    smoothing = 0.1

    temperature = 300
    V_readout = 0

    (edge_x, edge_y, charge), pot_descr = cc_3D.cc(n_eff=n_eff_0, bias=-20,
                                                   V_readout=V_readout,
                                                   temperature=temperature,
                                                   t_e_trapping=None,
                                                   t_h_trapping=None,
                                                   n_pixel_x=n_pixel_x,
                                                   n_pixel_y=n_pixel_y,
                                                   width_x=width_x,
                                                   width_y=width_y,
                                                   radius=radius, nD=nD,
                                                   resolution=resolution,
                                                   smoothing=smoothing)
 
    tools.save(obj=(edge_x, edge_y, charge),
               filename='charge_3D_unirrad')

    (edge_x, edge_y, charge), pot_descr = cc_3D.cc(n_eff=9.1e12, bias=-180,
                                                   V_readout=V_readout,
                                                   temperature=temperature,
                                                   t_e_trapping=2.6,
                                                   t_h_trapping=2.6,
                                                   n_pixel_x=n_pixel_x,
                                                   n_pixel_y=n_pixel_y,
                                                   width_x=width_x,
                                                   width_y=width_y,
                                                   radius=radius, nD=nD,
                                                   resolution=resolution,
                                                   smoothing=smoothing)

    tools.save(obj=(edge_x, edge_y, charge),
               filename='charge_3D_irrad')

    # Plot collected charge map
    _, _, charge_u = tools.load(filename='charge_3D_unirrad')
    plt.clf()
    plt.gca().set_aspect('equal')
    plt.gca().invert_yaxis()
    cmap = cm.get_cmap('inferno')
    cmap.set_bad('white')
    cmesh = plt.pcolormesh(edge_x, edge_y, charge/charge_u,
                           cmap=cmap, vmin=0, vmax=1.05)
    plt.title('Charge collection')
    plt.grid()
    cax = plt.gcf().add_axes(
        [plt.gca().get_position().xmax, 0.1, 0.05,
         plt.gca().get_position().ymax - plt.gca().get_position().ymin])
    plt.colorbar(cmesh, cax=cax, orientation='vertical')
    plt.grid()
    plt.savefig('CCE_3D.pdf', layout='tight')
    plt.show()
