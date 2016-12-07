''' Example that moves e-h pairs in a planar sensor.

    Calculates the induced charge from e-h pairs drifting
    through the silicon.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scarce import fields, plot, geometry, silicon, solver


def transient_planar():
    # Number of pixels influences how correct the field for the
    # center pixel(s) is due to more far away infinite boundary condition
    n_pixel = 9

    # Geometry of one pixel
    width = 50.
    thickness = 200.
    pitch = 45.

    n_eff = 5e12  # n_eff [cm^-3]
    temperature = 300

    # Potentials
    V_bias = -560.
    V_readout = 0.
    V_bi = -silicon.get_diffusion_potential(n_eff, temperature)

    # Create mesh of the sensor and stores the result
    # The created file can be viewed with any mesh viewer (e.g. gmsh)
    mesh = geometry.mesh_planar_sensor(
        n_pixel=n_pixel,
        width=width,
        thickness=thickness,
        resolution=300.,
        filename='planar_mesh_example.msh')

    # Numerically solve the laplace equation on the mesh
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

    w_potential = fields.calculate_planar_sensor_w_potential(
        mesh=mesh,
        width=width,
        pitch=pitch,
        n_pixel=n_pixel,
        thickness=thickness)

    min_x = float(mesh.getFaceCenters()[0, :].min())
    max_x = float(mesh.getFaceCenters()[0, :].max())
    pot_descr = fields.Description(potential,
                                   min_x=min_x,
                                   max_x=max_x,
                                   min_y=0,
                                   max_y=thickness)
    pot_w_descr = fields.Description(w_potential,
                                     min_x=min_x,
                                     max_x=max_x,
                                     min_y=0,
                                     max_y=thickness)

    # Start parameters of e-h pairs
    xx, yy = np.meshgrid(np.linspace(0, width, 2),  # x
                         np.linspace(thickness / 2., thickness / 2., 1),  # y
                         sparse=False)  # all combinations of x / y
    p0 = np.array([xx.ravel(), yy.ravel()])

    # Initial charge
    q0 = np.ones(p0.shape[1])

    # Time steps
    dt = 0.001  # [ns]
    n_steps = 3000
    t = np.arange(n_steps) * dt

    dd = solver.DriftDiffusionSolver(pot_descr, pot_w_descr, T=temperature)
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
    plot.get_planar_sensor_plot(fig=fig,
                                width=width,
                                pitch=pitch,
                                thickness=thickness,
                                n_pixel=n_pixel,
                                V_backplane=V_bias,
                                V_readout=V_readout,
                                potential_function=pot_descr.get_potential,
                                field_function=pot_descr.get_field,
                                depletion_function=pot_descr.get_depletion,
                                title='Planar sensor potential')

    # Create animation
    init, animate = plot.animate_drift_diffusion(
        fig, pe=traj_e, ph=traj_h, dt=dt)
    ani_time = 5.  # [s]
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
