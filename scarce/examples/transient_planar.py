''' Example that moves e-h pairs in a planar sensor.

    Calculates the induced charge from e-h pairs drifting
    through the silicon.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy import integrate

from scarce import plot, solver, sensor, tools


def transient_planar():
    # Sensor parameters
    n_eff = 1.45e12
    n_pixel = 9
    width = 50.
    pitch = 30.
    thickness = 200.
    smoothing = 0.05
    resolution = 251
    temperature = 300.
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
                                                  resolution=resolution,
                                                  # Might have to be adjusted
                                                  # when changing the geometry
                                                  smoothing=smoothing)

    # Start parameters of e-h pairs
    # Create 10 e-h pairs every 5 um in y
    xx, yy = np.meshgrid(np.linspace(0, width, 1),  # x
                         np.repeat(np.linspace(0.,
                                               thickness - 10,
                                               10), 1),
                         sparse=False)  # All combinations of x / y
    p0 = np.array([xx.ravel(), yy.ravel()])  # Position [um]

    # Initial charge set to 1
    q0 = np.ones(p0.shape[1])

    # Time steps
    dt = 0.001  # [ns]
    n_steps = 20000
    t = np.arange(n_steps) * dt

    dd = solver.DriftDiffusionSolver(pot_descr, pot_w_descr,
                                     T=temperature, diffusion=True,
                                     t_r=0.,
                                     save_frac=1)
    traj_e, traj_h, I_ind_e, I_ind_h, T, I_ind_tot, Q_ind_e_tot, Q_ind_h_tot = dd.solve(p0, q0, dt,
                                                              n_steps,
                                                              multicore=True)

    # Interpolate data to fixed time points for easier plotting
#     I_ind_e = tools.time_data_interpolate(T, I_ind_e, t, axis=0, fill_value=0.)
#     I_ind_h = tools.time_data_interpolate(T, I_ind_h, t, axis=0, fill_value=0.)
    I_ind_e[np.isnan(I_ind_e)] = 0.
    I_ind_h[np.isnan(I_ind_h)] = 0.
    Q_ind_e = integrate.cumtrapz(I_ind_e, T, axis=0, initial=0)
    Q_ind_h = integrate.cumtrapz(I_ind_h, T, axis=0, initial=0)
    Q_ind_t = integrate.cumtrapz(I_ind_tot, t, axis=0, initial=0)

    # Induced current of one e-h pair
#     plt.plot(T[:, 0], I_ind_e[:, 0],
#              '.-', color='blue', label='e', linewidth=2)
#     plt.plot(T[:, 0], I_ind_h[:, 0],
#              '.-', color='red', label='h', linewidth=2)
#     plt.plot(T[:, 0], I_ind_e[:, 0] + I_ind_h[:, 0], '.-', color='magenta',
#              label='e+h', linewidth=2)
#     plt.plot(t, I_ind_tot, '.-', color='black',
#              label='e+h', linewidth=2)
#     plt.ylabel('Induced current [a.u.]')
#     plt.grid()
#     plt.title('Signal of one drifting e-h pair, readout pixel')
#     plt.show()

    #     plt.plot(t, I_ind_tot,
#              color='black', label='Current all charges', linewidth=2)
#     plt.plot(t, np.cumsum(I_ind_tot) * dt, '--', linewidth=2,
#              color='black', label='all e+h')

    # Induced charge of one e-h pair
#     plt.plot(T[:, 0], Q_ind_e[:, 0], '.-',
#              color='blue', linewidth=2, label='e')
#     plt.plot(T[:, 0], Q_ind_h[:, 0], '.-',
#              color='red', linewidth=2, label='h')
    Q_ind_tot = Q_ind_e_tot + Q_ind_h_tot
    
    for i in range(10):
        plt.plot(T[:, i], Q_ind_e[:, i] + Q_ind_h[:, i],
                 '-', linewidth=1, label='e+h')
    
        plt.plot(plt.xlim(), [Q_ind_tot[i], Q_ind_tot[i]])
        print Q_ind_tot[i]

#     plt.plot(t, Q_ind_t, '.-', color='black',
#              label='e+h', linewidth=1)
    plt.legend(loc=0)
    plt.title('Induced charge of e-h pair, readout pixel')
    plt.xlabel('Time [ns]')
    plt.ylabel('Total induced charge [a.u.]')
    print Q_ind_e_tot, Q_ind_h_tot
    plt.show()
    
    

#     # Total induced charge of e-h pairs in readout pixel
#     plt.plot(t, Q_ind_e[:, ::2].sum(axis=1) / xx.shape[0], color='blue',
#              linewidth=2, label='e')
#     plt.plot(t, Q_ind_h[:, ::2].sum(axis=1) / xx.shape[0], color='red',
#              linewidth=2, label='h')
#     plt.plot(t, (Q_ind_e[:, ::2] + Q_ind_h[:, ::2]).sum(axis=1) / xx.shape[0],
#              color='magenta', linewidth=2, label='e+h')
#     plt.legend(loc=0)
#     plt.xlabel('Time [ns]')
#     plt.ylabel('Total induced charge [a.u.]')
#     plt.grid()
#     plt.title('Induced charge of MIP, readout pixel')
#     plt.show()
# 
#     # Mean induced charge of e-h pairs in readout pixel
#     plt.plot(t, Q_ind_e[:, 1::2].sum(axis=1) / xx.shape[0], color='blue',
#              linewidth=2, label='e')
#     plt.plot(t, Q_ind_h[:, 1::2].sum(axis=1) / xx.shape[0], color='red',
#              linewidth=2, label='h')
#     plt.plot(t, (Q_ind_e[:, 1::2] + Q_ind_h[:, 1::2]).sum(axis=1) / xx.shape[0],
#              color='magenta', linewidth=2, label='e+h')
#     plt.legend(loc=0)
#     plt.xlabel('Time [ns]')
#     plt.ylabel('Total induced charge [a.u.]')
#     plt.grid()
#     plt.title('Induced charge of MIP, neighbouring pixel')
#     plt.show()

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
