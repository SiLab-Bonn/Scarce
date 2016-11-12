import numpy as np

from scarce import silicon


def get_weighting_potential(x, y, D, S, W=None, is_planar=True):
    """ Planar sensor:
        From Nuclear Instruments and Methods in Physics Research A 535 (2004)
        554-557, with correction from wbar = pi*w/2/D to wbar = pi*w/D with:

        x [um] is the offset from the middle of the electrode
        y [um] the position in the sensor
        D [um] the sensor thickness
        S [um] the pixel pitch
        W [um] the electrode width

        3D sensor:
        Weighting potential for two cylinders with:
        D [um] distance between columns
        S [um] is the radius
        W number of readout columns
    """

    # Wheighting potential for one pixel
    if is_planar:
        if not W:  # Special case: 100% fill factor, amnalytic solution exists
            W = S

        xbar = np.pi * x / D
        ybar = np.pi * (y - D) / D
        wbar = np.pi * S / D
        return -1. / np.pi * (np.arctan(np.tan(ybar / 2) * np.tanh((xbar + wbar / 2.) / 2.)) -
                              np.arctan(np.tan(ybar / 2) * np.tanh((xbar - wbar / 2.) / 2.)))
    else:
        R = S
        D = D / 2.  # D is the total distance between the columns
        a = np.sqrt(D * D - R * R)
        Phi_w = 1. / (4 * np.arccosh(D / R)) * \
            np.log(((x - a) ** 2 + y ** 2) / ((x + a) ** 2 + y ** 2)) + 0.5

        # Stability
        Phi_w = np.ma.masked_where(
            np.sqrt((x + D) * (x + D) + y * y) < R, Phi_w)
        Phi_w = np.ma.masked_where(
            np.sqrt((x - D) * (x - D) + y * y) < R, Phi_w)
        Phi_w = np.ma.masked_where(Phi_w < 0., Phi_w)
        Phi_w = np.ma.masked_where(Phi_w > 1., Phi_w)

        return Phi_w


def get_weighting_field(x, y, D, S, is_planar=True):
    """ From Nuclear Instruments and Methods in Physics Research A 535 (2004)
        554-557, with correction from wbar = pi*w/2/D to wbar = pi*w/D
        with x [um] is the position in the sensor [0:thickness], y [um] the offset from the
        middle of the electrode, D [um] the sensor thickness and S [um] the
        eletrode width. The field is calculated from the drivation of the
        potential in x and y.
    """

    if is_planar:
        xbar = np.pi * x / D
        ybar = np.pi * (y - D) / D
        wbar = np.pi * S / D

        # Not easy to find a more simple form
        denom = (np.cosh(1. / 2. * (wbar - 2. * xbar)) + np.cos(ybar)) * \
            (np.cosh(1. / 2. * (wbar + 2. * xbar)) + np.cos(ybar)) * D

        E_x = - np.sin(ybar) * np.sinh(wbar / 2.) * np.sinh(xbar) / denom

        E_y = np.sinh(
            wbar / 2.) * (np.cosh(wbar / 2.) + np.cos(ybar) * np.cosh(xbar)) / denom

        return E_x, E_y
    else:
        # 3D sensor:
        # From the analytical derivation of the get_weighting_potential function
        # Weighting potential for two cylinders with:
        # S [um] is the radius
        # D [um] distance between columns

        R = S
        D = D / 2.
        a = np.sqrt(D * D - R * R)

        E_x = a / (np.arccosh(D / R)) * (a ** 2 - x ** 2 + y ** 2) / \
            (((a - x) ** 2 + y ** 2) * ((a + x) ** 2 + y ** 2))
        E_y = -2 * a / (np.arccosh(D / R)) * (x * y) / \
            (((a - x) ** 2 + y ** 2) * ((a + x) ** 2 + y ** 2))

        E_x = np.ma.masked_where(np.sqrt((x + D) * (x + D) + y * y) < R, E_x)
        E_x = np.ma.masked_where(np.sqrt((x - D) * (x - D) + y * y) < R, E_x)
        E_y = np.ma.masked_where(np.sqrt((x + D) * (x + D) + y * y) < R, E_y)
        E_y = np.ma.masked_where(np.sqrt((x - D) * (x - D) + y * y) < R, E_y)

        return -E_x, -E_y


def get_electric_field(x, y, V_bias, n_eff, D, S=None, is_planar=True):
    """ Calculates the 2D electric field E_x, E_y [V/um]

    Planar sensor:
        Calculates the field E_y[V/um], E_x = 0 in a planar sensor as a
        function of the position x between the electrodes [um],
        the bias voltage V_bias [V], the effective doping
        concentration n_eff [cm^-3] and the sensor Width D [um].
        The analytical function from the detector book p. 93 is used.

    3D sensor:
        Calculates the field E_x/E_y [V/um] in a 3d sensor as a function of the position
        x,y between the electrodes [um], the bias Voltage V_bias [V], the effective
        doping concentration n_eff [cm^-3], the electrode distance D [um] and radius R [um].
        So far the same field like the weighting field is used --> space charge is ignored.
    """

    if is_planar:
        if S:
            raise NotImplementedError(
                'The electrode width cannot be set, only full fill factor supported!')
        V_dep = silicon.get_depletion_voltage(n_eff, D)  # Depletion voltage
        a = (V_bias - V_dep) / D
        b = -2. * V_dep / (D ** 2)
        E_y = a - b * y
        return np.zeros_like(E_y), E_y
    else:
        E_x, E_y = get_weighting_field(x, y, D, S, is_planar=False)
        E_x = E_x * V_bias
        E_y = E_y * V_bias

        return E_x, E_y
