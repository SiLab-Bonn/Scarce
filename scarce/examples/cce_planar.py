''' Example that moves e-h pairs in a planar sensor.

    Calculates the induced charge from e-h pairs drifting
    through the silicon.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm
from scarce import plot, solver, geometry, silicon, fields, tools

import time
def timing(f):
    def wrap(*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        print '%s function took %0.3f ms' % (f.func_name, (time2-time1)*1000.0)
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
    # Create n_pairs e-h pairs every 5 um in x and y
    n_pairs = 1
    grid_x = 5  # grid spacing in x
    grid_y = 5  # grid spacing in y
    x_bins = int(width / grid_x)
    y_bins = int(thickness / grid_y)
    print x_bins, y_bins

    range_x=(-width / 2., width / 2.)
    range_y=(0, thickness)

    # Create e-h pairs in the pixel, avoid charge carriers on boundaries
    # e.g. x = -width / 2 or y = 0
    xx, yy = np.meshgrid(np.linspace(range_x[0] + grid_x / 2.,
                                     range_x[1] - grid_x / 2., x_bins),
                         np.repeat(np.linspace(range_y[0] + grid_y / 2.,
                                               range_y[1] - grid_y / 2., y_bins),
                                   n_pairs),  # 10 e-h per position
                         sparse=False)  # All combinations of x / y
    p0 = np.array([xx.ravel(), yy.ravel()])  # Position [um]
  
    # Initial charge set to 1
    q_start = 1
    q0 = np.ones(p0.shape[1]) * q_start
  
    # Time steps
    dt = 0.001  # [ns]
    n_steps = 20000

    dd = solver.DriftDiffusionSolver(pot_descr, pot_w_descr,
                                     T=temperature, diffusion=True, save_frac=1)
    traj_e, traj_h, I_ind_e, I_ind_h, t, _ = dd.solve(p0, q0, dt, n_steps,
                                                multicore=True)
 
    # Trajectory at t=0 is start position
    pos_0 = traj_e[0]

    # Last recorded integrated charge is total induced charge
    q_ind = np.cumsum(I_ind_e + I_ind_h, axis=0)[-1] * np.diff(t)[0]
    
#     q_ind = q0
#     pos_0 = p0
#     
#     print pos_0.shape
#     print q_ind.shape
#     
#     print pos_0[:, 0]
#     print q_ind[0]
    
    data = np.vstack((pos_0[0], pos_0[1], q_ind)).T
#     print data.shape
#     print data

    print 'q_ind<=0.1', data[q_ind<=0.1]
    
#     plt.clf()
#     plt.plot(t, np.cumsum(I_ind_e, axis=0)[:, 0] * dt, color='red',
#              label='Electrons')
#     plt.plot(t, np.cumsum(I_ind_h, axis=0)[:, 0] * dt, color='red',
#              label='Holes')
#     plt.plot(t, np.cumsum(I_ind_e + I_ind_h, axis=0)[:, 0] * dt,
#              color='magenta', label='Sum')
#     plt.show()

    n_bins_c = 100
    H, edges = np.histogramdd(sample=data,
                   bins=(x_bins, y_bins, n_bins_c),
                   range=((range_x[0], range_x[1]),
                          (range_y[0], range_y[1]),
                          (0., 2*q_start + 0.5 * q_start / n_bins_c)))
    
#     print range_x[0], range_x[1]
#     print range_y[0], range_y[1]

#     print (np.nansum(H, axis=2) == 0).shape
#     print q_ind.shape
#     raise
    
#     charge_pos = np.ma.masked_array(data=np.zeros(shape=(x_bins, y_bins)),mask=np.zeros(shape=(x_bins, y_bins)))
    charge_pos = np.zeros(shape=(x_bins, y_bins))
    sel = (np.sum(H, axis=2) != 0)

    weights = (edges[2][:-1] + edges[2][1:]) / 2.
    
    charge_pos[sel] = np.average(H, axis=2, weights=weights)[sel] * weights.sum() / np.sum(H, axis=2)[sel]
    
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
    cmesh = plt.pcolormesh(edges[0], edges[1], charge_pos.T, cmap=cmap, vmin=0, vmax=1.)#np.max(charge_pos))
    plt.grid()
    cax = plt.gcf().add_axes([plt.gca().get_position().xmax, 0.1, 0.05,
                              plt.gca().get_position().ymax - plt.gca().get_position().ymin])
    plt.colorbar(cmesh, cax=cax, orientation='vertical')
    plt.grid()
    plt.show()
#     
#     plt.clf()
#     plt.plot(t, np.cumsum(I_ind_e, axis=0)[:, 0] * dt, color='red',
#              label='Electrons')
#     plt.plot(t, np.cumsum(I_ind_h[:, 0], axis=0) * dt, color='red',
#              label='Holes')
#     plt.plot(t, np.cumsum(I_ind_e[:, 0] + I_ind_h[:, 0], axis=0) * dt,
#              color='magenta', label='Sum')
#     plt.show()
# 
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
#     
#     plot.get_planar_cce_plot(fig=fig,
#                                 width=width,
#                                 pitch=pitch,
#                                 thickness=thickness,
#                                 n_pixel=n_pixel,
#                                 
#                                 title='Planar sensor potential')
#   
#     print traj_e.shape
#     print q_ind.shape
#     print np.count_nonzero(q_ind<=0.01)

    # Create animation
    init, animate = plot.animate_drift_diffusion(fig,
                                                 pe=traj_e,
                                                 ph=traj_h, dt=dt)
    ani_time = 5.  # [ns]
    frames = 30
    ani = animation.FuncAnimation(fig=fig, func=animate,
                                  frames=np.arange(1, traj_h.shape[0],
                                                   traj_h.shape[0] / frames),
                                  interval=ani_time / frames * 1000.,
                                  blit=True, init_func=init,
                                  repeat_delay=ani_time / 5.)
  
    # ani.save('Example_planar_drift.gif', dpi=80, writer='imagemagick')
    plt.show()


if __name__ == '__main__':
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")
    transient_planar()
