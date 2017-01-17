''' Example that moves e-h pairs in a planar sensor.

    Calculates the induced charge from e-h pairs drifting
    through the silicon.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scarce import plot, solver, sensor


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

    # Create sensor
    pot_w_descr, pot_descr = sensor.planar_sensor(n_eff=n_eff,
                                                  V_bias=V_bias,
                                                  V_readout=V_readout,
                                                  temperature=temperature,
                                                  n_pixel=n_pixel,
                                                  width=width,
                                                  pitch=pitch,
                                                  thickness=thickness,
                                                  resolution=100.)

    # Start parameters of e-h pairs
    xx, yy = np.meshgrid(np.linspace(0, width, 2),  # x
                         np.linspace(3. * thickness / 4.,
                                     3. * thickness / 4., 2),
                         sparse=False)  # All combinations of x / y
    p0 = np.array([xx.ravel(), yy.ravel()])  # Position [um]

    # Initial charge set to 1
    q0 = np.ones(p0.shape[1])

    # Time steps
    dt = 0.001  # [ns]
    n_steps = 10000
    t = np.arange(n_steps) * dt

    dd = solver.DriftDiffusionSolver(pot_descr, pot_w_descr,
                                     T=temperature, diffusion=False)
    traj_e, traj_h, Q_ind_e, Q_ind_h, I_ind_e, I_ind_h = dd.solve(p0, q0, dt,
                                                                  n_steps)

    plt.plot(t, Q_ind_e[:, 0], color='blue', label='Electrons')
    plt.plot(t, Q_ind_h[:, 0], color='red', label='Holes')
    plt.plot(t, Q_ind_e[:, 0] + Q_ind_h[:, 0], color='magenta', label='Sum')
    plt.legend(loc=0)
    plt.xlabel('Time [ns]')
    plt.ylabel('Total induced charge [a.u.]')
    plt.twinx(plt.gca())
    plt.plot(t, I_ind_e[:, 0], '--', color='blue', label='Electrons')
    plt.plot(t, I_ind_h[:, 0], '--', color='red', label='Holes')
    plt.plot(t, I_ind_e[:, 0] + I_ind_h[:, 0], '--', color='magenta',
             label='Sum')
    plt.ylabel('Induced current [a.u.]')
    plt.grid()
    plt.title('Signal of drifting e-h pairs, readout pixel')
    plt.show()

    plt.plot(t, Q_ind_e[:, 1], color='blue', label='Electrons')
    plt.plot(t, Q_ind_h[:, 1], color='red', label='Holes')
    plt.plot(t, Q_ind_e[:, 1] + Q_ind_h[:, 1], color='magenta', label='Sum')
    plt.legend(loc=0)
    plt.xlabel('Time [ns]')
    plt.ylabel('Total induced charge [a.u.]')
    plt.twinx(plt.gca())
    plt.plot(t, I_ind_e[:, 1], '--', color='blue', label='Electrons')
    plt.plot(t, I_ind_h[:, 1], '--', color='red', label='Holes')
    plt.plot(t, I_ind_e[:, 1] + I_ind_h[:, 1], '--', color='magenta',
             label='Sum')
    plt.ylabel('Induced current [a.u.]')
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
