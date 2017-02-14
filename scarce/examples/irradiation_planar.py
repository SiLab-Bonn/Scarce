''' Example that shows the transient planar sensor signal after irradiation.
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

from scarce import silicon, solver, sensor, tools


def transient_irrad():
    # For CCE important parameters

    fluence = 5e15  # Neq/cm2
    V_bias = -1000.
    n_eff_0 = 1.7e12

    # Calculate effective doping concentration after irradiation
    n_eff = silicon.get_eff_acceptor_concentration(fluence=fluence,
                                                   n_eff_0=n_eff_0,
                                                   is_ntype=True,
                                                   is_oxygenated=True)

    # Planar sensor parameters
    width = 50.
    thickness = 200.
    temperature = 300.
    pitch = 45.
    n_pixel = 9
    V_readout = 0.
    resolution = 200

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
                                                  # when changing the
                                                  # geometry
                                                  smoothing=0.01)

    # Start parameters of e-h pairs
    # Create 10 e-h pairs every 5 um in y
    xx, yy = np.meshgrid(np.linspace(0, width, 1),  # x
                         np.repeat(np.linspace(0., thickness,
                                               thickness / 5.), 10),
                         sparse=False)  # All combinations of x / y
    p0 = np.array([xx.ravel(), yy.ravel()])  # Position [um]

    # Initial charge set to 1
    q0 = np.ones(p0.shape[1])

    # Time steps
    dt = 0.001  # [ns]
    n_steps = 3000
    t = np.arange(n_steps) * dt

    t_e_trapping = silicon.get_trapping(
        fluence=fluence, is_electron=True, paper=1)
    t_h_trapping = silicon.get_trapping(
        fluence=fluence, is_electron=False, paper=1)

    dd = solver.DriftDiffusionSolver(pot_descr, pot_w_descr,
                                     T=temperature, diffusion=True)
    dd_irr = solver.DriftDiffusionSolver(pot_descr, pot_w_descr,
                                         T=temperature, diffusion=True,
                                         t_e_trapping=t_e_trapping,
                                         t_h_trapping=t_h_trapping)
    _, _, I_ind_e, I_ind_h, T, _ = dd.solve(p0, q0, dt, n_steps)
    _, _, I_ind_e_irr, I_ind_h_irr, T_irr, _ = dd_irr.solve(p0, q0, dt,
                                                            n_steps)

    # Interpolate data to fixed time points for easier plotting
    I_ind_e = tools.time_data_interpolate(T, I_ind_e, t, axis=0, fill_value=0.)
    I_ind_h = tools.time_data_interpolate(T, I_ind_h, t, axis=0, fill_value=0.)
    I_ind_e[np.isnan(I_ind_e)] = 0.
    I_ind_h[np.isnan(I_ind_h)] = 0.
    I_ind_e_irr = tools.time_data_interpolate(
        T_irr, I_ind_e_irr, t, axis=0, fill_value=0.)
    I_ind_h_irr = tools.time_data_interpolate(
        T_irr, I_ind_h_irr, t, axis=0, fill_value=0.)
    I_ind_e_irr[np.isnan(I_ind_e_irr)] = 0.
    I_ind_h_irr[np.isnan(I_ind_h_irr)] = 0.
    Q_ind_e = integrate.cumtrapz(I_ind_e, t, axis=0, initial=0)
    Q_ind_h = integrate.cumtrapz(I_ind_h, t, axis=0, initial=0)
    Q_ind_e_irr = integrate.cumtrapz(I_ind_e_irr, t, axis=0, initial=0)
    Q_ind_h_irr = integrate.cumtrapz(I_ind_h_irr, t, axis=0, initial=0)
    plt.plot(t, Q_ind_e.sum(axis=1) / xx.shape[0], color='blue',
             label='Electrons, depl.')
    plt.plot(t, Q_ind_h.sum(axis=1) / xx.shape[0], color='red',
             label='Holes, depl.')
    plt.plot(t, (Q_ind_e.sum(axis=1) +
                 Q_ind_h.sum(axis=1)) / xx.shape[0], color='magenta',
             label='Sum, depl.')
    plt.plot(t, Q_ind_e_irr.sum(axis=1) / xx.shape[0], '--', color='blue',
             label='Electrons, depl. + trap.')
    plt.plot(t, Q_ind_h_irr.sum(axis=1) / xx.shape[0], '--', color='red',
             label='Holes, depl. + trap.')
    plt.plot(t, (Q_ind_e_irr.sum(axis=1) +
                 Q_ind_h_irr.sum(axis=1)) / xx.shape[0], '--', color='magenta',
             label='Sum, depl. + trap.')
    plt.legend(loc=0)
    plt.xlabel('Time [ns]')
    plt.ylabel('Total induced charge [a.u.]')
    plt.grid()
    plt.title('Induced charge of MIP in planar sensor, readout pixel')
    plt.show()


if __name__ == '__main__':
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")
    transient_irrad()
