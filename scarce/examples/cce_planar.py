''' Example that moves e-h pairs in a planar sensor.

    Calculates the induced charge from e-h pairs drifting
    through the silicon.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm
from scipy import integrate
from scarce import plot, solver, geometry, silicon, fields, tools

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


def transient_planar():
    # Sensor parameters
    n_eff = 1.7e12
    width = 50.
    thickness = 200.
    temperature = 300.
    pitch = 45.
    n_pixel = 9
    V_bias = -80.
    V_readout = 0.
    resolution = 300
    smoothing = 0.05

    # Create mesh of the sensor and stores the result
    # The created file can be viewed with any mesh viewer (e.g. gmsh)
    mesh = geometry.mesh_planar_sensor(
        n_pixel=n_pixel,
        width=width,
        thickness=thickness,
        resolution=resolution)

    min_x = float(mesh.getFaceCenters()[0, :].min())
    max_x = float(mesh.getFaceCenters()[0, :].max())

    # Set um resolution grid
    nx = width * n_pixel
    ny = thickness

    V_bi = -silicon.get_diffusion_potential(n_eff, temperature)
    # Numerically solve the Laplace equation on the mesh
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
    # Numerically solve the Poisson equation on the mesh
    w_potential = fields.calculate_planar_sensor_w_potential(
        mesh=mesh,
        width=width,
        pitch=pitch,
        n_pixel=n_pixel,
        thickness=thickness)
    pot_w_descr = fields.Description(w_potential,
                                     min_x=min_x,
                                     max_x=max_x,
                                     min_y=0,
                                     max_y=thickness,
                                     nx=nx,
                                     ny=ny,
                                     smoothing=smoothing)

    # Start parameters of e-h pairs
    # Create n_pairs e-h pairs
    n_pairs = 200
    grid_x = 10  # grid spacing in x in um
    grid_y = 10  # grid spacing in y in um
    x_bins = int(width / grid_x)
    y_bins = int(thickness / grid_y)

    print x_bins, y_bins

    range_x = (-width / 2., width / 2.)
    range_y = (0, thickness)

    # Create e-h pairs in the pixel, avoid charge carriers on boundaries
    # e.g. x = -width / 2 or y = 0
    xx, yy = np.meshgrid(np.linspace(range_x[0] + grid_x / 2.,
                                     range_x[1] - grid_x / 2., x_bins),
                         np.repeat(np.linspace(range_y[0] + grid_y / 2.,
                                               range_y[1] - grid_y / 2., y_bins),
                                   n_pairs),  # 10 e-h per position
                         sparse=False)  # All combinations of x / y

#     xx, yy = np.meshgrid(np.linspace(17.5,
#                                      17.5, 1),
#                          np.repeat(np.linspace(137.5,
#                                                137.5, 1),
#                                    n_pairs),  # 10 e-h per position
#                          sparse=False)  # All combinations of x / y

    p0 = np.array([xx.ravel(), yy.ravel()])  # Position [um]

    # Initial charge set to 1
    q_start = 1.
    q0 = np.ones(p0.shape[1]) * q_start
    q_max = q_start * 1.05

    # Time steps
    dt = 0.001  # [ns]
    n_steps = 25000

    t = np.linspace(0, n_steps * dt, 1000)

    dd = solver.DriftDiffusionSolver(pot_descr, pot_w_descr,
                                     T=temperature, diffusion=True, save_frac=50)
    traj_e, traj_h, I_ind_e, I_ind_h, T, _, Q_ind_e_tot, Q_ind_h_tot = dd.solve(p0, q0, dt, n_steps,
                                                      multicore=True)

    # Trajectory at t=0 is start position
    pos_0 = traj_e[0]

    # Interpolate data to fixed time points for easier plotting
#     I_ind_e = tools.time_data_interpolate(T, I_ind_e, t, axis=0, fill_value=0.)
#     I_ind_h = tools.time_data_interpolate(T, I_ind_h, t, axis=0, fill_value=0.)
    I_ind_e[np.isnan(I_ind_e)] = 0.
    I_ind_h[np.isnan(I_ind_h)] = 0.
    Q_ind_e = integrate.cumtrapz(I_ind_e, T, axis=0, initial=0)
    Q_ind_h = integrate.cumtrapz(I_ind_h, T, axis=0, initial=0)

    # Last index with data (time != nan)
    index = np.nanargmax(T, axis=0)
    y = np.indices(index.shape)
    # Last recorded integrated charge is total induced charge
    q_ind = Q_ind_e[index, y][0] + Q_ind_h[index, y][0]
    q_ind = Q_ind_e_tot + Q_ind_h_tot

    # Histogram charge per start position
    data = np.vstack((pos_0[0], pos_0[1], q_ind)).T

    n_bins_c = 100
    H, edges = np.histogramdd(sample=data,
                              bins=(x_bins, y_bins, n_bins_c),
                              range=((range_x[0], range_x[1]),
                                     (range_y[0], range_y[1]),
                                     (0., q_max)))

#     sel = q_ind > 1.03
    # Induced charge of one e-h pair
#     plt.plot(T, Q_ind_e, '-',
#              color='blue', linewidth=2, label='e')
#     plt.plot(T, Q_ind_h, '-',
#              color='red', linewidth=2, label='h')
#     plt.plot(T, Q_ind_e + Q_ind_h,
#              '-', color='magenta', linewidth=2, label='e+h')
#     plt.legend(loc=0)
#     plt.xlabel('Time [ns]')
#     plt.ylabel('Total induced charge [a.u.]')

#     plt.twinx()
#     plt.plot(t[:, sel], Q_ind_e[:, sel], '.-',
#              color='blue', linewidth=2, label='e')
#     plt.plot(t[:, sel], Q_ind_h[:, sel], '.-',
#              color='red', linewidth=2, label='h')
#     plt.plot(t[:, sel], Q_ind_e[:, sel] + Q_ind_h[:, sel],
#              '.-', color='magenta', linewidth=2, label='e+h')

#     plt.show()
#     print pos_0[0][sel], pos_0[1][sel]

#     print range_x[0], range_x[1]
#     print range_y[0], range_y[1]

#     print (np.nansum(H, axis=2) == 0).shape
#     print q_ind.shape
#     raise

#     charge_pos = np.ma.masked_array(data=np.zeros(shape=(x_bins, y_bins)),mask=np.zeros(shape=(x_bins, y_bins)))
    charge_pos = np.zeros(shape=(x_bins, y_bins))
    sel = (np.sum(H, axis=2) != 0)

    weights = (edges[2][:-1] + edges[2][1:]) / 2.

    charge_pos[sel] = np.average(
        H, axis=2, weights=weights)[sel] * weights.sum() / np.sum(H, axis=2)[sel]

#     print np.average(H, axis=2, weights=weights)
#     raise
    # charge_pos = np.nansum(H, axis=2)
#     charge_pos.mask[~sel] = 1

    print edges[0].shape
    print edges[1].shape
    print 'sel.shape', sel.shape
    print 'charge_pos.shape', charge_pos.shape
    print np.max(charge_pos)
    print '((edges[0][:-1] + edges[0][1:]) / 2.)', ((edges[0][:-1] + edges[0][1:]) / 2.)

    print ((edges[0][:-1] + edges[0][1:]) / 2.)[sel[0][0]]
    print ((edges[1][:-1] + edges[1][1:]) / 2.)[sel[0][1]]

    print np.count_nonzero(np.average(H, axis=2, weights=weights)), data.shape[0]
#     raise

#     print 'sel', sel
#     raise
    # X, Y = np.meshgrid((edges[0][:-1] + edges[0][1:])/2., (edges[1][:-1] + edges[1][1:])/2.)

    plt.clf()
    plt.gca().set_aspect('equal')
    plt.gca().invert_yaxis()

    cmap = cm.get_cmap('inferno')
    cmap.set_bad('white')
    # np.max(charge_pos))
    cmesh = plt.pcolormesh(
        edges[0], edges[1], charge_pos.T, cmap=cmap, vmin=0, vmax=q_max)
    plt.grid()
    cax = plt.gcf().add_axes([plt.gca().get_position().xmax, 0.1, 0.05,
                              plt.gca().get_position().ymax - plt.gca().get_position().ymin])
    plt.colorbar(cmesh, cax=cax, orientation='vertical')
    plt.grid()
    plt.show()

    # Plot numerical result in 2D with particle animation
    fig = plt.figure()
    plot.get_planar_sensor_plot(fig=fig,
                                width=width,
                                pitch=pitch,
                                thickness=thickness,
                                n_pixel=n_pixel,
                                V_backplane=V_bias,
                                V_readout=V_readout,
                                pot_func=pot_descr.get_potential,
                                field_func=pot_descr.get_field,
                                depl_func=pot_descr.get_depletion,
                                title='Planar sensor potential')

    # Create animation
    frames = 50
    init, animate = plot.animate_drift_diffusion(fig, T=T, pe=traj_e,
                                                 ph=traj_h, dt=t.max() /
                                                 frames,
                                                 n_steps=frames)
    ani = animation.FuncAnimation(fig=fig, func=animate,
                                  blit=True, init_func=init, frames=frames,
                                  interval=5000 / frames)

    # ani.save('Example_planar_drift.gif', dpi=80, writer='imagemagick')
    plt.show()


if __name__ == '__main__':
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")
    transient_planar()
