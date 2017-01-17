''' Example that shows the transient signal after irradiation.

    Planar and 3D sensor signal is simulated for different fluences.
'''

import numpy as np
import matplotlib.pyplot as plt
from scarce import silicon, solver, sensor


def transient_irrad():
    # Planar sensor parameters
    n_eff = 1.7e12
    width = 50.
    thickness = 200.
    temperature = 300.
    pitch = 45.
    n_pixel = 9
    V_bias = -1000.
    V_readout = 0.

    # Create sensor
    planar_pot_w, planar_pot = sensor.planar_sensor(n_eff=n_eff,
                                             V_bias=V_bias,
                                             V_readout=V_readout,
                                             temperature=temperature,
                                             n_pixel=n_pixel,
                                             width=width,
                                             pitch=pitch,
                                             thickness=thickness,
                                             resolution=200.)

    # Start parameters of e-h pairs
    # Create a e-h pairs every 5 um in y
    xx, yy = np.meshgrid(np.linspace(0, width, 1),  # x
                         np.linspace(0, thickness, thickness/5.),
                         sparse=False)  # All combinations of x / y
    p0 = np.array([xx.ravel(), yy.ravel()])  # Position [um]

    # Initial charge set to 1
    q0 = np.ones(p0.shape[1])

    # Time steps
    dt = 0.001  # [ns]
    n_steps = 3000
    t = np.arange(n_steps) * dt

    t_e_trapping = silicon.get_trapping(
        fluence=1e15, is_electron=True, paper=1)
    t_h_trapping = silicon.get_trapping(
        fluence=1e15, is_electron=False, paper=1)

    dd = solver.DriftDiffusionSolver(planar_pot, planar_pot_w,
                                     T=temperature, diffusion=False)
    dd_irrad = solver.DriftDiffusionSolver(planar_pot, planar_pot_w,
                                           T=temperature, diffusion=False,
                                           t_e_trapping=t_e_trapping,
                                           t_h_trapping=t_h_trapping)
    _, _, Q_ind_e, Q_ind_h, _, _ = dd.solve(p0, q0, dt, n_steps)

    _, _, Q_ind_e_irr, Q_ind_h_irr, _, _ = dd_irrad.solve(p0, q0, dt, n_steps)

    plt.plot(t, Q_ind_e.sum(axis=1) / xx.shape[0], color='blue',
             label='Electrons')
    plt.plot(t, Q_ind_h.sum(axis=1) / xx.shape[0], color='red',
             label='Holes')
    plt.plot(t, (Q_ind_e.sum(axis=1) +
                 Q_ind_h.sum(axis=1)) / xx.shape[0], color='magenta',
             label='Sum')
    plt.plot(t, Q_ind_e_irr.sum(axis=1) / xx.shape[0], '--', color='blue',
             label='Electrons irrad.')
    plt.plot(t, Q_ind_h_irr.sum(axis=1) / xx.shape[0], '--', color='red',
             label='Holes irrad.')
    plt.plot(t, (Q_ind_e_irr.sum(axis=1) +
                 Q_ind_h_irr.sum(axis=1)) / xx.shape[0], '--', color='magenta',
             label='Sum irrad.')
    plt.legend(loc=0)
    plt.xlabel('Time [ns]')
    plt.ylabel('Total induced charge [a.u.]')
    plt.grid()
    plt.title('Induced charge of MIP, readout pixel')
    plt.show()


if __name__ == '__main__':
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")
    transient_irrad()
