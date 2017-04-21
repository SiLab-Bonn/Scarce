''' Example that moves e-h pairs in a 3D sensor.

    Calculates the induced charge from e-h pairs drifting
    through the silicon.

    .. WARNING::
       The 3D field is low between pixels. Thus diffusion
       should be activated to leave this field minima quickly
       to give reasonable results.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy import integrate

from scarce import plot, solver, sensor


def transient_3D():
    n_pixel_x, n_pixel_y = 3, 3
    width_x = 250.
    width_y = 50.
    radius = 6.
    nD = 2  # Number of columns per pixel
    n_eff = 1e12
    T = 300
    V_bias = -20.
    V_readout = 0.

    pot_w_descr, pot_descr, geom_descr = sensor.sensor_3D(n_eff=n_eff,
                                                          V_bias=V_bias,
                                                          V_readout=V_readout,
                                                          temperature=T,
                                                          n_pixel_x=n_pixel_x,
                                                          n_pixel_y=n_pixel_y,
                                                          width_x=width_x,
                                                          width_y=width_y,
                                                          radius=radius,
                                                          nD=nD,
                                                          resolution=80,
                                                          smoothing=0.1)

    # Start parameters of e-h pairs
    # 10 pairs per position
    xx, yy = np.meshgrid(np.linspace(-width_x / 2., width_x / 2., 25),  # x
                         np.repeat(np.linspace(-width_y / 2.,
                                               width_y / 2., 5), 10),  # y
                         sparse=False)  # All combinations of x / y
    p0 = np.array([xx.ravel(), yy.ravel()])

    # Initial charge
    q0 = np.ones(p0.shape[1])

    # Time steps
    dt = 0.001  # [ns]
    n_steps = 10000
    t = np.arange(n_steps) * dt

    dd = solver.DriftDiffusionSolver(pot_descr, pot_w_descr,
                                     geom_descr=geom_descr,
                                     T=T,
                                     diffusion=True)
    traj_e, traj_h, I_ind_e, I_ind_h, T, I_ind_tot, Q_ind_e_tot, Q_ind_h_tot = dd.solve(p0, q0, dt, n_steps)

    I_ind_e[np.isnan(I_ind_e)] = 0.
    I_ind_h[np.isnan(I_ind_h)] = 0.
    Q_ind_e = integrate.cumtrapz(I_ind_e, T, axis=0, initial=0)
    Q_ind_h = integrate.cumtrapz(I_ind_h, T, axis=0, initial=0)

    for i in (5, 13):
        plt.plot(T[:, i], Q_ind_e[:, i], color='blue',
                 label='Electrons')
        plt.plot(T[:, i], Q_ind_h[:, i], color='red',
                 label='Holes')
        plt.plot(T[:, i], (Q_ind_e + Q_ind_h)[:, i],
                 color='magenta', label='Sum')
        plt.plot(plt.xlim(), [(Q_ind_e_tot + Q_ind_h_tot)[i]] * 2, '-', color='magenta')
        plt.legend(loc=0)
        plt.xlabel('Time [ns]')
        plt.ylabel('Charge normalized to 1')
        plt.grid()
        plt.title('Induced charge of drifting e-h pairs, start pos. %d/%d um' %
                  (p0[0, i], p0[1, i]))
        plt.show()

    # Plot numerical result in 2D with particle animation
    fig = plt.figure()
    plot.get_3D_sensor_plot(fig, width_x, width_y,
                            radius, nD,
                            n_pixel_x, n_pixel_y,
                            V_bias=V_bias, V_readout=V_readout,
                            pot_func=pot_descr.get_potential_smooth,
                            field_func=pot_descr.get_field,
                            # Comment in if you want to see the mesh
                            mesh=None,  # potential.mesh,
                            title='Potential and field of a 3D sensor, '\
                            '%dx%d pixel matrix, numerical solution' % \
                            (n_pixel_x, n_pixel_y))

    # Create animation
    frames = 50
    init, animate = plot.animate_drift_diffusion(fig, T=T, pe=traj_e,
                                                 ph=traj_h, dt=t.max() /
                                                 frames,
                                                 n_steps=frames)
    ani = animation.FuncAnimation(fig=fig, func=animate,
                                  blit=True, init_func=init, frames=frames,
                                  interval=5000 / frames)

    # ani.save('Example_3D_drift.gif', dpi=80, writer='imagemagick')
    plt.show()


if __name__ == '__main__':
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")
    transient_3D()
