""" Collection of (semi-) empiric functions describing silicon properties.
"""

import numpy as np
from scipy import constants

from scarce import constant


def get_depletion_depth(V_bias, n_eff, temperature):
    ''' Depletion depth [um] of silcon with an effective doping concentration n_eff [10^12 /cm^3]
        at a temperature [K] and reverse bias V_bias [V].

        Check/citetation of formulars needed!'''

    # Relative permitivity of silicon
    epsilon_r = constant.epsilon_s / constants.epsilon_0

    V_bi = get_diffusion_potential(n_eff=n_eff, temperature=temperature)

    # Depletion depth [m]
    return np.sqrt(
        2. * constants.epsilon_0 * epsilon_r * (V_bi + V_bias) / (constants.elementary_charge * n_eff * 10 ** 6))


def get_depletion_voltage(n_eff, distance):
    """ This function returns the full depletion Voltage [V] as a function
        of the effective doping concentration Neff [10^12 /cm^-3] and the
        distance between electrodes [um]. Formular is standard and for example
        used in:
        G. Kramberger et al., Nucl. Inst. And. Meth. A 476 (2002) 645-651
        'Determination of effective trapping times for electrons and holes
        in irradiated silicon'.
    """

    return constants.elementary_charge * n_eff / (constant.epsilon_s) * distance ** 2. / 2 * 1e6


def get_diffusion_potential(n_eff, temperature):
    ''' Diffusion potential [V] as a function of the effective doping
        concentration n_eff [10^12 / cm^3] and the temperature [K].

        Check/citetation of formulars needed!'''

    # [10^12/cm^3] empirical fit at K = 300 range
    N_i = 9.38 * 10 ** 7. * \
        (temperature / 300.) ** 2 * np.exp(-6884. / temperature)

    # Diffusion potential in thermal equilibrium [V]
    return constants.Boltzmann * temperature / constants.elementary_charge * np.log(n_eff ** 2. / N_i ** 2)


def get_eff_acceptor_concentration(fluence, n_eff_0, is_ntype, is_oxygenated):
    """ Calculate the effective acceptor concentration [10^12 cm^-3] of irradiated n and
        p-type silicon with and without oxygen enriched as a function of the
        fluence [10^12 cm^-2]. The data can be desribed by different line fits.
        The parameters were extracted from a plot taken for n-type silicon from:
        CERN/LHCC 2000-009 LEB Status Report/RD48 31 Dec. 1999
        and for p-type silicon from:
        RD50 Status Report 2006 CERN-LHCC-2007-005, p. 4-6

        Due to the difference in the data for different technologies
        a rather large error on the propotionality factor of 10% is assumed.
    """

    fluence = np.atleast_1d(fluence)

    # Constants
    if is_ntype:  # n-type silicon
        # [cm^-1], slope before type inversion
        beta_n = -5.5e-2
        # [cm^-1], slope after type inversion, standart silicon FZ
        beta_n_fz = 2.83e-2
        # [cm^-1], slope after type inversion, oxygenated silicon FZ
        beta_n_ofz = 0.94e-2
        Phi_inv = -n_eff_0 / beta_n  # Inversion point N_eq [10^14 cm^-2]

        # After type inversion
        if is_oxygenated:  # Oxygenated silicon
            n_eff = beta_n_ofz * (fluence - Phi_inv)
        else:  # Standart silicon
            n_eff = beta_n_fz * (fluence - Phi_inv)

        # Before type inversion
        n_eff[fluence < Phi_inv] = n_eff_0 + \
            beta_n * fluence[fluence < Phi_inv]

    else:  # p-type silicons
        beta_n_fz = 2.3e-2  # [cm^-1]
        beta_n_ofz = 2.7e-2  # [cm^-1]
        if is_oxygenated:  # Oxygenated silicon
            n_eff = n_eff_0 + beta_n_ofz * fluence
        else:  # Standart silicon
            n_eff = n_eff_0 + beta_n_fz * fluence

    return n_eff


def get_free_path(fluence, e_field, temperature, is_electron=True):

    # Calculate the mean free path in cm of the charge carriers from the trapping
    # propability and the velocity. The trapping propability is a
    # function of the fluence and the velocity is a function
    # of the electric field and the temperature. The electric field itself
    # depends on the electrode geometry and the bias voltage.
    velocity = get_mobility(e_field, temperature, is_electron) * e_field
    trapping_time = get_trapping(fluence, is_electron, paper=1)

    return velocity * trapping_time * 10e-9


def get_mobility(e_field, temperature, is_electron):
    ''' Calculate the mobility [cm^2/Vs] of charge carriers in silicon from
        the electrical field (E [V/cm]) the temperature (T [K]) and the charge
        carrier type (isElectron [0/1] otherwise hole).

        Formular derived from measured data of high purity silicon and the
        corresponding fit function parameters are used here.
        From:
        C. Jacononi et al., Solid state electronics, 1977, vol 20., p. 87
        'A review of some charge transport properties of silicon'
        Note: the doping concentration is irrelevant for < 10^16/cm^3
    '''

    if is_electron:
        v_m = 1.53e9 * temperature ** (-0.87)  # [cm/s]
        E_c = 1.01 * temperature ** 1.55  # [V/cm]
        beta = 2.57e-2 * temperature ** 0.66
    else:  # Holes
        v_m = 1.62e8 * temperature ** (-0.52)  # [cm/s]
        E_c = 1.24 * temperature ** 1.68  # [V/cm]
        beta = 0.46 * temperature ** 0.17

    mu = v_m / E_c / (1. + (np.abs(e_field) / E_c) ** beta) ** (1. / beta)

    return mu


def get_resistivity(n_eff, is_n_type=True, temperature=300, e_field=1e3):
    ''' Calculate the resitivity from:
        The effective doping concentration n_eff [10^12 / cm^3]
        the mobility [cm^2/Vs]
        for n- and p-type silicon.

        The mobility istself is a
        function of the temperature [K] and the electric field [V/cm].
        From http://ecee.colorado.edu/~bart/book/mobility.htm

        TODO: If you take the mobility[E_field] equation seriously, then there is no constant
        resitivity since the mobility depends also on the electric field. For low E-Fields <= 1000 V/cm
        the mobility is independent of the E flied and thus the resistivity. Likely this parameter
        is always given in low field approximation?! Source needed!
    '''

    mobility = get_mobility(e_field, temperature, is_electron=is_n_type)

    return 1. / (constants.e * n_eff * mobility * 1e12)


def get_trapping(fluence, is_electron, paper=1):
    ''' Calculate the trapping time tr (e^-(tr) in ns) of charge carriers in
        silicon as a function of the fluence. There was also a dependence on
        the temperature measured, that is omitted here!
    '''

    if paper == 1:
        # Newer paper where the finding was that there is no difference
        # in electron / hole trapping propability. Electron trapping is similar
        # to Kramberger et. al. but hole trapping much lower.
        # IEEE TRANSACTIONS ON NUCLEAR SCIENCE, VOL. 51, NO. 6, DECEMBER 2004
        # 'Measurement of Trapping Time Constants in
        # Proton-Irradiated Silicon Pad Detectors'

        if is_electron:
            beta = 5.13e-16  # [cm^2/ns]
            # beta_error = 0.16e-16  # [cm^2/ns]
        else:
            beta = 5.04e-16  # [cm^2/ns]
            # beta_error = 0.18e-16  # [cm^2/ns]

        tr = 1. / (fluence * beta)
        # tr_error = 1. / (fluence * beta ** 2) * beta_error

        return tr

    elif paper == 2:
        # Oldest most cited reference, with irradiation to 2 e 14 only.
        # Calculates the trapping time tr (e^-(tr) in ns) of charge carriers in
        # silicon at a temperature of T = 263 K as a function of the
        # fluence (in 10^12 Neq/cm^2). This was measured with a laser and planar
        # silicon sensors with a fluence up to 2*10^14 Neq/cm^2. There was a
        # linear behaviour between fluence and the effective trapping
        # propability measured intepended of the silicon type (oxygenated or
        # not and with different doping concentrations) from:
        # G. Kramberger et al., Nucl. Inst. And. Meth. A 476 (2002) 645-651
        # 'Determination of effective trapping times for electrons and holes
        # in irradiated silicon'

        if is_electron:
            beta = 4.2e-16  # [cm^2/ns]
            # beta_error = 0.3e-16  # [cm^2/ns]
        else:
            beta = 6.1e-16  # [cm^2/ns]
            # beta_error = 0.3e-16  # [cm^2/ns]

        tr = 1. / (fluence * beta)
        # tr_error = 1. / (fluence * beta ** 2) * beta_error

        return tr
    else:
        raise RuntimeError('Unknown paper selected!')
