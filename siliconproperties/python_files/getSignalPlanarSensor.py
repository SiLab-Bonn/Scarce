import numpy as np
from scipy import constants

from siliconproperties.python_files.getEffAcceptorConcentration import get_eff_acceptor_concentration
from siliconproperties.python_files.getTrapping import get_trapping
from siliconproperties.python_files.getElectricField import get_field
from siliconproperties.python_files.getWeightingPotential import get_weighting_potential
from siliconproperties.python_files.getWeightingField import get_weighting_field
from siliconproperties.python_files.getMobility import get_mobility


def get_signal_planar_sensor(x_track, q_in, D, S, n_eff_0, is_ntype, is_oxygenated, V_bias, fluence, temperatur, t, dt, N):
    ''' Function to return the signal of a planar silicon sensor as a function of
        the time. The parameters are:
        x_track: offset from the electrode mean [um]
        q_in: deposited charge [e-h pairs]
        D: sensor thickness [um]
        S: electrode width [um]
        n_eff_0: doping concentration without fluence [10^12 Neq/cm^3]
        isNtype: = 1 for ntype sensor
        isOxygenated: = 1 for oxygenated sensor
        V_bias: bias of the sensor [V]
        fluence: radiation dose [10^12 Neq/cm^3]
        temperatur: temerature [K]
        t: time points of simulation [ns]
        dt: time step of simulation [ns] CHOOSE SMALL ENOUGH!
        N: number of quasi particles along the track
    '''

    # Constants
    Q_e = constants.e  # Electron charge [C]
    e_h = q_in  # deposited charge in electron hole pairs
    x_0 = x_track  # offset from the middle of the sensor pixel (x-direction) [um]

    n_eff = get_eff_acceptor_concentration(fluence, n_eff_0, is_ntype, is_oxygenated)

    if fluence:
        tr_e = get_trapping(fluence, is_electron=True, paper=1)  # [ns]
        tr_h = get_trapping(fluence, is_electron=False, paper=1)  # [ns]

    e_h = e_h / N

    # Start positions of the electrons
    if N == 1:
        y_e = np.atleast_1d(D / 2.)
    else:
        y_e = np.linspace(0, D, N)

    y_h = y_e.copy()  # start positions of the electron holes
    x_e = np.ones_like(y_h) * x_0  # Positions of the electrons in x
    x_h = x_e.copy()  # Positions of the electron holes in x
    Q_ind_e = np.zeros_like(y_e)  # start induced charge from electrons
    Q_ind_h = np.zeros_like(y_e)  # start induced charge from holes

    Q_ind_e_vec = np.zeros(shape=(y_e.shape[0], t.shape[0]))
    Q_ind_h_vec = np.zeros(shape=(y_e.shape[0], t.shape[0]))
    y_e_vec = np.zeros(shape=(y_e.shape[0], t.shape[0]))
    y_h_vec = np.zeros(shape=(y_e.shape[0], t.shape[0]))
    v_e_vec = np.zeros(shape=(y_e.shape[0], t.shape[0]))
    v_h_vec = np.zeros(shape=(y_e.shape[0], t.shape[0]))
    E_q_y_e_vec = np.zeros(shape=(y_e.shape[0], t.shape[0]))
    E_q_y_h_vec = np.zeros(shape=(y_e.shape[0], t.shape[0]))
    E_w_y_e_vec = np.zeros(shape=(y_e.shape[0], t.shape[0]))
    E_w_y_h_vec = np.zeros(shape=(y_e.shape[0], t.shape[0]))
    Phi_w_y_e_vec = np.zeros(shape=(y_e.shape[0], t.shape[0]))
    Phi_w_y_h_vec = np.zeros(shape=(y_e.shape[0], t.shape[0]))

    for index, i in enumerate(t):  # time loop
        # Electric, weighting field a charge carrier positions
        # Electrons:
        _, E_q_y_e = get_field(x=x_e,  # Electric field [V/um]
                               y=y_e,
                               V_bias=V_bias,
                               n_eff=n_eff,
                               D=D,
                               S=None,  # FIXME
                               is_planar=True)

        Phi_w_e = get_weighting_potential(x=x_e,
                                          y=y_e,
                                          D=D,
                                          S=S,
                                          is_planar=True)

        _, E_w_y_e = get_weighting_field(x=x_e,
                                         y=y_e,
                                         D=D,
                                         S=S,
                                         is_planar=True)

        # Holes
        _, E_q_y_h = get_field(x=x_h,  # Electric field [V/um]
                               y=y_h,
                               V_bias=V_bias,
                               n_eff=n_eff,
                               D=D,
                               S=None,  # FIXME
                               is_planar=True)

        Phi_w_h = get_weighting_potential(x=x_h,
                                          y=y_h,
                                          D=D,
                                          S=S,
                                          is_planar=True)

        _, E_w_y_h = get_weighting_field(x=x_h,
                                         y=y_h,
                                         D=D,
                                         S=S,
                                         is_planar=True)

        # Movement:

        # Electrons
        v_e = - E_q_y_e * get_mobility(E_q_y_e * 1e5,  # Velocity [um/ns]
                                       temperature=temperatur,
                                       is_electron=True)
        dy_e = v_e * dt
        y_e = y_e + dy_e
        # Boundaries
        E_w_y_e[y_e > D] = 0
        E_q_y_e[y_e > D] = 0
        v_e[y_e >= D] = 0
        y_e[y_e > D] = D

        # Holes
        v_h = E_q_y_h * get_mobility(E_q_y_h * 1e5,  # Velocity [um/ns]
                                     temperature=temperatur,
                                     is_electron=False)
        dy_h = v_h * dt
        y_h = y_h + dy_h
        # Boundaries
        E_w_y_h[y_h < 0] = 0
        E_q_y_h[y_h < 0] = 0
        v_h[y_h <= 0] = 0
        y_h[y_h < 0] = 0

        # Induced charge calculation with trapping
        # electrons
        dQ_e = np.zeros_like(y_e)
        dQ_e[y_e < D] = e_h * E_w_y_e[y_e < D] * dy_e[y_e < D]
        if fluence:
            dQ_e[y_e < D] *= np.exp(-i / tr_e)

        # Holes
        dQ_h = np.zeros_like(y_h)
        dQ_h[y_h > 0] = -e_h * E_w_y_h[y_h > 0] * dy_h[y_h > 0]
        if fluence:
            dQ_h[y_h > 0] *= np.exp(-i / tr_h)

        # Sum up
        Q_ind_e += dQ_e
        Q_ind_h += dQ_h

        # Data for plotting
        Q_ind_e_vec[:, index] = Q_ind_e
        Q_ind_h_vec[:, index] = Q_ind_h
        y_e_vec[:, index] = y_e
        v_e_vec[:, index] = v_e
        y_h_vec[:, index] = y_h
        v_h_vec[:, index] = v_h
        E_q_y_e_vec[:, index] = E_q_y_e
        E_q_y_h_vec[:, index] = E_q_y_h
        E_w_y_e_vec[:, index] = E_w_y_e
        E_w_y_h_vec[:, index] = E_w_y_h
        Phi_w_y_e_vec[:, index] = Phi_w_e
        Phi_w_y_h_vec[:, index] = Phi_w_h

    print Q_ind_e_vec
    # It is already summed?
    # Q_tot = sum(Q_ind_e_vec, 1) + sum(Q_ind_h_vec, 1)


if __name__ == '__main__':
    import matplotlib.pylab as plt

    dt = 1e-3
    t = np.linspace(0., 6., 6. / dt)

    get_signal_planar_sensor(x_track=0,  # Offset from the middle of the sensor pixel (x-direction) [um]
                             q_in=23e3,  # Deposited charge in electron hole pairs
                             D=300,  # Detector width [um]
                             S=100,  # Pixel width [um]
                             n_eff_0=1.,  # Doping concentration [10^12 /cm^3]
                             is_ntype=True,
                             is_oxygenated=True,
                             V_bias=80,
                             fluence=0,
                             temperatur=248,
                             t=t,
                             dt=dt,
                             N=1)
