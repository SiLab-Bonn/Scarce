import numpy as np
import progressbar
from scipy import constants

from siliconproperties.python_files.getEffAcceptorConcentration import get_eff_acceptor_concentration
from siliconproperties.python_files.getTrapping import get_trapping
from siliconproperties.python_files.getElectricField import get_field
from siliconproperties.python_files.getWeightingPotential import get_weighting_potential
from siliconproperties.python_files.getWeightingField import get_weighting_field
from siliconproperties.python_files.getMobility import get_mobility
from siliconproperties.python_files.getDepletionVoltage import get_depletion_voltage


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
        tr_e = get_trapping(fluence * 1e12, is_electron=True, paper=2)  # [ns]
        tr_h = get_trapping(fluence * 1e12, is_electron=False, paper=2)  # [ns]

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

    progress_bar = progressbar.ProgressBar(widgets=['', progressbar.Percentage(), ' ', progressbar.Bar(marker='*', left='|', right='|'), ' ', progressbar.AdaptiveETA()], maxval=len(t), term_width=80)
    progress_bar.start()

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
        e_in_boundary = np.logical_and(y_e <= D, y_e >= 0)
        E_w_y_e[~e_in_boundary] = 0
        E_q_y_e[~e_in_boundary] = 0
        v_e[~e_in_boundary] = 0
        y_e[~e_in_boundary] = D
        E_w_y_e[~e_in_boundary] = 0
        E_q_y_e[~e_in_boundary] = 0
        v_e[~e_in_boundary] = 0
        y_e[~e_in_boundary] = 0

        # Holes
        v_h = E_q_y_h * get_mobility(E_q_y_h * 1e5,  # Velocity [um/ns]
                                     temperature=temperatur,
                                     is_electron=False)
        dy_h = v_h * dt
        y_h = y_h + dy_h
        # Boundaries
        h_in_boundary = np.logical_and(y_h <= D, y_h >= 0)
        E_w_y_h[~h_in_boundary] = 0
        E_q_y_h[~h_in_boundary] = 0
        v_h[~h_in_boundary] = 0
        y_h[~h_in_boundary] = 0
        E_w_y_h[~h_in_boundary] = 0
        E_q_y_h[~h_in_boundary] = 0
        v_h[~h_in_boundary] = 0
        y_h[~h_in_boundary] = D

        # Induced charge calculation with trapping
        # electrons
        dQ_e = np.zeros_like(y_e)
        dQ_e[e_in_boundary] = e_h * E_w_y_e[e_in_boundary] * dy_e[e_in_boundary]
        if fluence:
            dQ_e[e_in_boundary] *= np.exp(-i / tr_e)

        # Holes
        dQ_h = np.zeros_like(y_h)
        dQ_h[h_in_boundary] = -e_h * E_w_y_h[h_in_boundary] * dy_h[h_in_boundary]
        if fluence:
            dQ_h[h_in_boundary] *= np.exp(-i / tr_h)

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
        
        progress_bar.update(index)
    progress_bar.finish()

    # It is already summed?
    Q_tot_e = Q_ind_e_vec.sum(axis=0)
    Q_tot_h = Q_ind_h_vec.sum(axis=0)
    Q_tot = Q_tot_e + Q_tot_h
    
    return y_e_vec, y_h_vec, Q_tot, Q_tot_e, Q_tot_h


if __name__ == '__main__':
    import matplotlib.pylab as plt

    dt = 1e-3
    t = np.linspace(0., 6., 6. / dt)
    
    CCE = []
    
    q_in = 23e3
    thickness = 200
    fluences = np.arange(10) * 1000.
    n_eff_0 = 1.7
    is_ntype=True
    is_oxygenated=True
    
#     for fluence in fluences:
#         n_eff = get_eff_acceptor_concentration(fluence, n_eff_0, is_ntype, is_oxygenated)
#         bias = get_depletion_voltage(n_eff, thickness)[0]
#         print fluence / 1000., bias
#     raise

#     4.6 e 11 - 1.4 e 12

    n_eff = get_eff_acceptor_concentration(5000, n_eff_0=0.5, is_ntype=False, is_oxygenated=True)
    print get_depletion_voltage(n_eff, 65)[0]
    
    raise

    n_eff = get_eff_acceptor_concentration(1000., n_eff_0, is_ntype, is_oxygenated)
    bias = get_depletion_voltage(n_eff, thickness)[0]
    voltages = [bias * 1.1] #range(int(bias), 1000, 100)
    CCE = []
 
    for voltage in voltages:
        y_e_vec, y_h_vec, Q_tot, Q_tot_e, Q_tot_h = get_signal_planar_sensor(x_track=0,  # Offset from the middle of the sensor pixel (x-direction) [um]
                                         q_in=q_in,  # Deposited charge in electron hole pairs
                                         D=thickness,  # Detector width [um]
                                         S=50,  # Pixel width [um]
                                         n_eff_0=n_eff_0,  # Doping concentration [10^12 /cm^3]
                                         is_ntype=is_ntype,
                                         is_oxygenated=is_oxygenated,
                                         V_bias=voltage,
                                         fluence=1000,  # 1e15
                                         temperatur=248,
                                         t=t,
                                         dt=dt,
                                         N=20)
#         print np.max(-Q_tot) / q_in
        CCE.append(np.max(-Q_tot) / q_in)
         
    print voltages, CCE
    plt.plot(voltages, CCE)
    plt.show()

#     for fluence in fluences:
#         n_eff = get_eff_acceptor_concentration(fluence, n_eff_0, is_ntype, is_oxygenated)
#         bias = get_depletion_voltage(n_eff, thickness)[0] * 1.1
# 
#         y_e_vec, y_h_vec, Q_tot, Q_tot_e, Q_tot_h = get_signal_planar_sensor(x_track=0,  # Offset from the middle of the sensor pixel (x-direction) [um]
#                                          q_in=q_in,  # Deposited charge in electron hole pairs
#                                          D=thickness,  # Detector width [um]
#                                          S=50,  # Pixel width [um]
#                                          n_eff_0=n_eff_0,  # Doping concentration [10^12 /cm^3]
#                                          is_ntype=is_ntype,
#                                          is_oxygenated=is_oxygenated,
#                                          V_bias=bias,
#                                          fluence=fluence,
#                                          temperatur=248,
#                                          t=t,
#                                          dt=dt,
#                                          N=20)
        

#         plt.plot(t, -Q_tot, label='Q_tot')
#         plt.plot(t, -Q_tot_e, label='Q_tot_e')
#         plt.plot(t, -Q_tot_h, label='Q_tot_h')
#         plt.legend(loc=0)
#         plt.show()
        
#         CCE.append(np.max(-Q_tot) / q_in)
    
#     plt.plot(fluences, CCE)
#     plt.plot(t, y_e_vec[0, :], label='electrons')
#     plt.plot(t, y_h_vec[0, :], label='holes')

#     plt.plot(t, b)
#     plt.legend(loc=0)
#     plt.show()
