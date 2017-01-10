''' Example that moves e-h pairs in a planar sensor.

    Calculates the induced charge from e-h pairs drifting
    through the silicon.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scarce import (fields, plot, geometry,
                    silicon, solver)


def transient_3D():
    # Number of pixels influences how correct the field for the
    # center pixel(s) is due to more correct boundary condition
    n_pixel_x, n_pixel_y = 3, 3

    # Geometry of one pixel
    width_x = 250.
    width_y = 50.
    radius = 6.
    nD = 2  # Number of columns per pixel

    n_eff = 1e12  # n_eff [cm^-3]
    temperature = 300

    # Potentials
    V_bias = -20.
    V_readout = 0.
    V_bi = -silicon.get_diffusion_potential(n_eff, temperature)

    resolution = 150.

    mesh = geometry.mesh_3D_sensor(width_x=width_x,
                                   width_y=width_y,
                                   n_pixel_x=n_pixel_x,
                                   n_pixel_y=n_pixel_y,
                                   radius=radius,
                                   nD=nD,
                                   resolution=resolution)

    potential = fields.calculate_3D_sensor_potential(mesh,
                                                     width_x,
                                                     width_y,
                                                     n_pixel_x,
                                                     n_pixel_y,
                                                     radius,
                                                     nD,
                                                     n_eff,
                                                     V_bias,
                                                     V_readout,
                                                     V_bi)

    w_potential = fields.calculate_3D_sensor_w_potential(mesh,
                                                         width_x,
                                                         width_y,
                                                         n_pixel_x,
                                                         n_pixel_y,
                                                         radius,
                                                         nD=nD)

    # Describe the 3D sensor array
    geom_descr = geometry.SensorDescription3D(width_x, width_y,
                                              n_pixel_x, n_pixel_y,
                                              radius, nD)
    min_x, max_x, min_y, max_y = geom_descr.get_array_corners()

    # Describe the result to be able to obtain field/potential at any point in
    # space
    pot_descr = fields.Description(potential,
                                   min_x=min_x,
                                   max_x=max_x,
                                   min_y=min_y,
                                   max_y=max_y,
                                   nx=width_x * n_pixel_x * 10.,  # um res.
                                   ny=width_y * n_pixel_y * 10.,  # um res.
                                   smoothing=0.1
                                   )
    pot_w_descr = fields.Description(w_potential,
                                     min_x=min_x,
                                     max_x=max_x,
                                     min_y=min_y,
                                     max_y=max_y,
                                     nx=width_x * n_pixel_x * 10.,
                                     ny=width_y * n_pixel_y * 10.,
                                     smoothing=0.1
                                     )

    # Start parameters of e-h pairs
    xx, yy = np.meshgrid(np.linspace(-width_x / 3., width_x / 3., 4),  # x
                         np.linspace(-width_y / 3., width_y / 3., 4),  # y
                         sparse=False)  # all combinations of x / y
    p0 = np.array([xx.ravel(), yy.ravel()])

    # Initial charge
    q0 = np.ones(p0.shape[1])

    # Time steps
    dt = 0.0001  # [ns]
    n_steps = 30000
    t = np.arange(n_steps) * dt

    dd = solver.DriftDiffusionSolver(pot_descr, pot_w_descr,
                                     T=temperature, geom_descr=geom_descr)
    traj_e, traj_h, Q_ind_e, Q_ind_h, Q_ind_tot = dd.solve(p0, q0, dt, n_steps)

    plt.plot(t, Q_ind_e[:, 0], label='Electrons')
    plt.plot(t, Q_ind_h[:, 0], label='Holes')
    plt.plot(t, Q_ind_tot[:, 0], label='Sum')
    plt.legend(loc=0)
    plt.xlabel('Time [ns]')
    plt.ylabel('Charge normalized to 1')
    plt.grid()
    plt.title('Induced charge of drifting e-h pairs, readout pixel')
    plt.show()

    plt.clf()
    plt.plot(t, Q_ind_e[:, 1], label='Electrons')
    plt.plot(t, Q_ind_h[:, 1], label='Holes')
    plt.plot(t, Q_ind_tot[:, 1], label='Sum')
    plt.legend(loc=0)
    plt.xlabel('Time [ns]')
    plt.ylabel('Charge normalized to 1')
    plt.grid()
    plt.title('Induced charge of drifting e-h pairs, neighbouring pixel')
    plt.show()

    # Plot numerical result in 2D with particle animation
    fig = plt.figure()
    plot.get_3D_sensor_plot(fig, width_x, width_y,
                            radius, nD,
                            n_pixel_x, n_pixel_y,
                            V_bias=V_bias + V_bi, V_readout=V_readout,
                            potential_function=pot_descr.get_potential_smooth,
                            field_function=pot_descr.get_field,
                            # Comment in if you want to see the mesh
                            mesh=None,  # potential.mesh,
                            title='Potential and field of a 3D sensor, '\
                            '%dx%d pixel matrix, numerical solution' % \
                            (n_pixel_x, n_pixel_y))

    # Create animation
    init, animate = plot.animate_drift_diffusion(
        fig, pe=traj_e, ph=traj_h, dt=dt)
    ani_time = 5.  # [s]
    frames = 100
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
    transient_3D()
