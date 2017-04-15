''' Example that moves e-h pairs in a planar sensor.

    Calculates the induced charge from e-h pairs drifting
    through the silicon.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm
from scipy import integrate
from scarce import plot, solver, geometry, silicon, fields, tools, analysis

import time
from matplotlib.pyplot import twinx


def timing(f):
    def wrap(*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        print '%s function took %0.3f ms' % (f.func_name, (time2 - time1) * 1000.0)
        return ret
    return wrap


def cc_efficiency(fluences, biases, temperatures, n_effs, t_e_trappings, t_h_trappings):
    n_pixel = 9
    width = 50.
    pitch = 30.
    thickness = 250.

    smoothing = 0.05
    resolution = 287

    try:
        n_effs = float(n_effs)
        n_effs = [n_effs] * len(n_effs)
    except TypeError:
        pass

    try:
        biases = float(biases)
        biases = [biases] * len(fluences)
    except TypeError:
        pass

    try:
        temperatures = float(temperatures)
        temperatures = [temperatures] * len(fluences)
    except TypeError:
        pass

    # Charge collection efficiency maps
    ccs = []

    # Create mesh of the sensor and stores the result
    mesh = geometry.mesh_planar_sensor(
        n_pixel=n_pixel,
        width=width,
        thickness=thickness,
        resolution=resolution)
    # Numerically solve the Poisson equation on the mesh
    # To create the weighting potential
    w_potential = fields.calculate_planar_sensor_w_potential(
        mesh=mesh,
        width=width,
        pitch=pitch,
        n_pixel=n_pixel,
        thickness=thickness)
    min_x = float(mesh.getFaceCenters()[0, :].min())
    max_x = float(mesh.getFaceCenters()[0, :].max())
    nx = width * n_pixel
    ny = thickness
    pot_w_descr = fields.Description(w_potential,
                                     min_x=min_x,
                                     max_x=max_x,
                                     min_y=0,
                                     max_y=thickness,
                                     nx=nx,
                                     ny=ny,
                                     smoothing=smoothing)

    pot_descr = None
    for i, fluence in enumerate(fluences):

        print 'fluence', fluence
        print 'n_eff', n_effs[i]
        print 'V_bias', biases[i]
        print 't_r', t_e_trappings[i], t_h_trappings[i]

        # Reset potenital description if recalc of field is needed
        # Is needed when n_eff changes
        if i > 0 and n_effs[i] != n_effs[i - 1]:
            pot_descr = None

        edge_x, edge_y, charge, pot_descr = charge_collected(n_eff=n_effs[i], t_e_trapping=t_e_trappings[i], t_h_trapping=t_h_trappings[i],
                                                             width=width, pitch=pitch, thickness=thickness, n_pixel=n_pixel,
                                                             temperature=temperatures[
                                                                 i], V_bias=biases[i],
                                                             pot_w_descr=pot_w_descr, mesh=mesh, smoothing=smoothing, pot_descr=pot_descr)
        # Plot collected charge map
        plt.clf()
        plt.gca().set_aspect('equal')
        plt.gca().invert_yaxis()
        cmap = cm.get_cmap('inferno')
        cmap.set_bad('white')
        cmesh = plt.pcolormesh(edge_x, edge_y, charge,
                               cmap=cmap, vmin=0, vmax=1.05)
        plt.grid()
        cax = plt.gcf().add_axes([plt.gca().get_position().xmax, 0.1, 0.05,
                                  plt.gca().get_position().ymax - plt.gca().get_position().ymin])
        plt.colorbar(cmesh, cax=cax, orientation='vertical')
        plt.grid()
        plt.title('Charge collection, fluence %1.2f neq_cm2' % fluence)
        ccs.append(charge)
        plt.savefig('CC_%d_%d.pdf' % (fluence, biases[i]), layout='tight')
        # plt.show()

        if i > 0:
            # Plot charge collection efficiency map
            plt.clf()
            plt.gca().set_aspect('equal')
            plt.gca().invert_yaxis()
            cmap = cm.get_cmap('inferno')
            cmap.set_bad('white')
            cmesh = plt.pcolormesh(edge_x, edge_y, charge / ccs[0] * 100.,
                                   cmap=cmap, vmin=0, vmax=100.)
            plt.grid()
            cax = plt.gcf().add_axes([plt.gca().get_position().xmax, 0.1, 0.05,
                                      plt.gca().get_position().ymax - plt.gca().get_position().ymin])
            plt.colorbar(cmesh, cax=cax, orientation='vertical')
            plt.grid()
            plt.title('Charge collection efficiency, '
                      'fluence %1.2f neq_cm2' % fluence)
            print 'CCE', fluence, charge.mean() / ccs[0].mean() * 100.
            plt.savefig('CCE_%d_%d.pdf' % (fluence, biases[i]), layout='tight')
            # plt.show()

    return ccs


def charge_collected(n_eff, t_e_trapping, t_h_trapping, width, pitch, thickness, n_pixel, temperature,
                     V_bias, pot_w_descr, mesh, smoothing, pot_descr=None):
    V_readout = 0.

    min_x = float(mesh.getFaceCenters()[0, :].min())
    max_x = float(mesh.getFaceCenters()[0, :].max())

    # Set um resolution grid
    nx = width * n_pixel
    ny = thickness

    V_bi = -silicon.get_diffusion_potential(n_eff, temperature)
    # Numerically solve the Laplace equation on the mesh
    if not pot_descr:
        potential = fields.calculate_planar_sensor_potential(
            mesh=mesh,
            width=width,
            pitch=pitch,
            n_pixel=n_pixel,
            thickness=thickness,
            n_eff=n_eff,
            V_bias=V_bias,
            V_readout=V_readout,
            V_bi=V_bi)
        pot_descr = fields.Description(potential,
                                       min_x=min_x,
                                       max_x=max_x,
                                       min_y=0,
                                       max_y=thickness,
                                       nx=nx,
                                       ny=ny,
                                       smoothing=smoothing)

    return analysis.get_charge_planar(width, thickness, pot_descr, pot_w_descr,
                                      t_e_trapping=t_e_trapping, t_h_trapping=t_h_trapping,
                                      grid_x=5, grid_y=5, n_pairs=10, dt=0.001,
                                      n_steps=20000, temperature=300), pot_descr

    #         # Plot numerical result in 2D with particle animation
    #         fig = plt.figure()
    #         plot.get_planar_sensor_plot(fig=fig,
    #                                     width=width,
    #                                     pitch=pitch,
    #                                     thickness=thickness,
    #                                     n_pixel=n_pixel,
    #                                     V_backplane=V_bias,
    #                                     V_readout=V_readout,
    #                                     pot_func=pot_descr.get_potential,
    #                                     field_func=pot_descr.get_field,
    #                                     depl_func=pot_descr.get_depletion,
    #                                     title='Planar sensor potential')
    #
    #         # Create animation
    #         frames = 50
    #         init, animate = plot.animate_drift_diffusion(fig, T=T, pe=traj_e,
    #                                                      ph=traj_h, dt=t.max() /
    #                                                      frames,
    #                                                      n_steps=frames)
    #         ani = animation.FuncAnimation(fig=fig, func=animate,
    #                                       blit=True, init_func=init, frames=frames,
    #                                       interval=5000 / frames)
    #
    #         # ani.save('Example_planar_drift.gif', dpi=80, writer='imagemagick')
    #         plt.show()


if __name__ == '__main__':
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")

    fluences = [1000.] * 4
    biases = [-200, -300, -400, -600]
    n_eff_0s = [1.475e12] * 4

    n_effs = []
    t_e_trappings = []
    t_h_trappings = []
    
    with open("data.txt", "a+") as myfile:
        myfile.append('Fluence\tBias\tN_eff_0\tN_eff\tt_e\tt_h\n')
    
        for i, fluence in enumerate(fluences):
            n_effs.append(silicon.get_eff_acceptor_concentration(fluence, n_eff_0s[i] / 1e12,
                                                                 is_ntype=True,
                                                                 is_oxygenated=True)[0] * 1e12)
            t_e_trappings.append(silicon.get_trapping(fluence * 1e12, is_electron=True,
                                                      paper=1) * (i + 1))
            t_h_trappings.append(silicon.get_trapping(fluence * 1e12, is_electron=False,
                                                      paper=1) * (i + 1))
    
            ccs = cc_efficiency(fluences=fluences,
                                biases=biases,
                                temperatures=300,
                                n_effs=n_effs,
                                t_e_trappings=t_e_trappings,
                                t_h_trappings=t_h_trappings)
            
        for i, fluence in enumerate(fluences):
            myfile.append('%d\t%d\t%1.2e\t%1.2e\t%1.2e\t%1.2e\n', (fluence,
                                                                   biases[i],
                                                                   n_eff_0s[i]))


    print biases
    print ccs_mean
#     ccs_mean = [100.0, 95.819278274882421, 91.492897261571542, 84.44507498549774, 78.023854736701296, 72.517359275661676, 67.581338385657943, 49.639626527411856, 38.39041632424324, 30.536375588272836, 24.835251542831479, 17.418036145128809, 13.384192805371708]
#     print biases
#     print ccs_mean
#     biases =  np.abs([-100, -200, -300, -400, -500, -600, -700, -800])
#     ccs_mean =  [30.315111552256575, 41.117291542671161, 48.451824125239519, 54.58876529456731, 58.059736778438165, 59.735573362868422, 60.508863282827484, 61.292900450661357]
    plt.clf()
    plt.plot(biases, ccs_mean, '.-')

    plt.show()
