""" Collection of (semi-) empiric functions describing particle interaction with silicon.
"""

import numpy as np
import os
from scipy import interpolate

from scarce import constant


def get_attenuation(energy, process='total'):
    """ Interpolates the data from http://physics.nist.gov/cgi-bin/Xcom
    to return the attenuation length for the given photon energy [MeV]
    in silicon for:
        - compton effect
        - photo effect
        - pair production
        - total
    """

    data = np.loadtxt(fname=os.path.join(constant.DATA_FOLDER,
                                         'AttenuationData.txt'),
                      dtype=np.float,
                      skiprows=11).T

    x = data[0]

    if 'compton' in process.lower():
        y = data[2]
    elif 'photo' in process.lower():
        y = data[3]
    elif 'pair' in process.lower():
        y = data[4] + data[5]
    elif 'total' in process.lower():
        y = data[6]
    else:
        raise AttributeError(
            'Unknown process %s. Try compton, photo, pair, total.', process)

    y *= constant.density  # cm2 / g -> 1 / cm

    # , kind='quadratic')
    interpol = interpolate.interp1d(
        x, y, bounds_error=False, fill_value=np.NaN)
    y_interpolated = interpol(energy)

    return y_interpolated


def get_dedx(betagamma, thickness, density_correction=True, restricted=True, MP=False):
    # FIXME STILL HAS BUGS; DOES NOT GIVE CORRECT VALUES!
    ''' Thickness [um] '''

    if restricted and not density_correction:
        raise NotImplementedError(
            'Restricted energy loss can only be calculated with density correction on')

    # Some derivations of the input to simplify the formulars
    betagamma2 = betagamma * betagamma
    gamma = (betagamma2 + 1) ** 0.5  # Calculate gamma from betagamma squared
    # Calculate beta from betagamma squared
    beta = (1 - 1. / (1 + betagamma2)) ** 0.5
    beta2 = beta * beta
    assert np.allclose(beta, (1 - 1. / (gamma * gamma)) ** 0.5)  # Cross check

    z = 1  # Charge of incident particle [e]
    mp = 938.27231              # mass of a proton [MeV/c**2]
    re = 2.817940325 * 10 ** -13     # electron radius [cm]
    Z = 14                      # atomic number [#nucleons]
    A = 28.0855                 # atomic mass of silicon absorber [amu = g/mol]
    Na = 6.0221415 * 10 ** 23        # avogrado number
    me = 0.510998918            # mass of an electron [MeV/c**2]
    I = 173 * 10 ** (-6)             # mean excitation energy [MeV]
    K = 4 * np.pi * Na * re ** 2 * me         # constant(?) [MeV]
    Alpha = 1 / 137.03599911      # fine structure constant

    # Density correction parameters for silicon, from
    # Atomic data and nuclear tables 30, 261-271 (1984)
    Omega_p = 31.06             # plasma energy in silicon [eV]
    C_d = -4.4351               # density correction parameter
    A_d = 0.14921               # density correction parameter
    M_d = 3.2546                # density correction parameter
    X_0_d = 0.2014              # density correction parameter
    X_1_d = 2.8715              # density correction parameter

    # Restricted Eloss, thin absorbers correction parameters
    # [MeV] for 200 um silicon from Nucl. Phys. B288 (1987) 681-716
    T_cut = 19e-3

    # MP Eloss, Landau-Vavilov Bichsel
    detector_density = 2.3290   # density of detector material [g cm**-3]
    # material length [g/cm**2]
    mat_length = thickness * 1e-4 * detector_density
    j = 0.198  # parameter [no unit]

    # Calculate dEdx
    # Maximum energy transfer possible in single collision
    Tmax = (2.0 * me * betagamma2) / \
        (1.0 + 2.0 * gamma * (me / mp) + (me / mp) ** 2)

    # dedx from semiempiric Bethe Bloch formular
    if not density_correction and not restricted:
        dedx = (-K * z ** 2 * Z) / (A * beta2) * \
            (0.5 * np.log((2 * me * betagamma2 * Tmax) / I ** 2) - beta2)
        return dedx

    # dEdx with density effect (with silicon parameters)
    # Density effect correction, from
    # PDG 30. Passage of particles through matter
    X = np.log(betagamma)  # Derivation of the input to simplify the formular
    # Calculate delta with three cases
    # Else case: zero for X < X_0
    delta = np.zeros_like(betagamma)
    # Case X >= X1: 2*ln 10 * X - C
    delta[X >= X_1_d] = 2 * np.log(10) * X[X >= X_1_d] + C_d
    selection = np.logical_and(X >= X_0_d, X < X_1_d)
    # Case X>=X_0_d & X<X_1_d: 2*log(10)*X + C_d + A_d*(X_1_d - X)**M_d
    delta[selection] = 2 * \
        np.log(10) * X[selection] + C_d + A_d * (X_1_d - X[selection]) ** M_d

    if density_correction and not restricted:
        dedx = (-K * z ** 2 * Z) / (A * beta2) * (0.5 *
                                                  np.log((2 * me * betagamma2 * Tmax) / I ** 2) - beta2 - delta / 2)
        return dedx

    # dedx with restricted Eloss correction
    if density_correction and restricted and not MP:
        dedx = (-K * z ** 2 * Z) / (A * beta2) * (0.5 * np.log((2 * me * betagamma2 * T_cut) / I ** 2) - beta2 /
                                                  2 * (1 + T_cut / Tmax) - delta / 2)  # dEdx with density effect and restriction on measured dEdx in thin absorbers
        return dedx

    # dedx MP Eloss correction
    if density_correction and restricted and MP:
        kappa = K / 2 * Z / A * mat_length / beta2
        dedx = kappa * (np.log(2 * me * betagamma2 / I) +
                        np.log(kappa / I) + j - beta2 - delta) / (mat_length)
        return -dedx
