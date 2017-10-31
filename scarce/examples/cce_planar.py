''' Example that moves e-h pairs in a planar sensor.

    Calculates the induced charge from e-h pairs drifting
    through the silicon.
'''

import matplotlib.pyplot as plt
from matplotlib import cm

from scarce import tools
from scarce.examples import cc_planar


if __name__ == '__main__':
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")

    n_pixel = 9
    width = 50.
    pitch = 30.
    thickness = 200.
    V_readout = 0.
    temperature = 300

    smoothing = 0.05
    resolution = 251

    # Unirradiated charge map
    (edge_x, edge_y, charge), pot_descr = cc_planar.cc(n_eff=1.475e12,
                                                       bias=-80,
                                                       V_readout=V_readout,
                                                       temperature=temperature,
                                                       t_e_trapping=None,
                                                       t_h_trapping=None,
                                                       n_pixel=n_pixel,
                                                       width=width,
                                                       pitch=pitch,
                                                       thickness=thickness,
                                                       resolution=resolution,
                                                       smoothing=smoothing)
    tools.save(obj=(edge_x, edge_y, charge),
               filename='charge_planar_unirrad')

    # Irradiated charge map
    (edge_x, edge_y, charge), pot_descr = cc_planar.cc(n_eff=8.8e12,
                                                       bias=-150,
                                                       V_readout=V_readout,
                                                       temperature=temperature,
                                                       t_e_trapping=4.82,
                                                       t_h_trapping=4.82,
                                                       n_pixel=n_pixel,
                                                       width=width,
                                                       pitch=pitch,
                                                       thickness=thickness,
                                                       resolution=resolution,
                                                       smoothing=smoothing)
    tools.save(obj=(edge_x, edge_y, charge),
               filename='charge_planar_irrad')

    # Plot CCE map
    _, _, charge_u = tools.load(filename='charge_planar_unirrad')
    plt.clf()
    plt.gca().set_aspect('equal')
    plt.gca().invert_yaxis()
    cmap = cm.get_cmap('inferno')
    cmap.set_bad('white')
    cmesh = plt.pcolormesh(edge_x, edge_y, charge / charge_u,
                           cmap=cmap, vmin=0, vmax=1.05)
    plt.title('Charge collection efficiency')
    plt.grid()
    cax = plt.gcf().add_axes(
        [plt.gca().get_position().xmax, 0.1, 0.05,
         plt.gca().get_position().ymax - plt.gca().get_position().ymin])
    plt.colorbar(cmesh, cax=cax, orientation='vertical')
    plt.grid()
    plt.savefig('CCE_planar_irrad.pdf', layout='tight')
    plt.show()
