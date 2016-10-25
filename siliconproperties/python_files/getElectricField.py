import numpy as np

from siliconproperties.python_files.getDepletionVoltage import get_depletion_voltage
from siliconproperties.python_files.getWeightingField import get_weighting_field


def get_field(x, y, V_bias, n_eff, D, S=None, is_planar=True):
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
            raise NotImplementedError('The electrode width cannot be set, only full fill factor supported!')
        V_dep = get_depletion_voltage(n_eff, D)  # Depletion voltage
        a = (V_bias - V_dep) / D
        b = -2. * V_dep / (D ** 2)
        E_y = a - b * y
        return np.zeros_like(E_y), E_y
    else:
        E_x, E_y = get_weighting_field(x, y, D, S, is_planar=False)
        E_x = E_x * V_bias
        E_y = E_y * V_bias

        return E_x, E_y


if __name__ == '__main__':
    import matplotlib.pylab as plt
    from siliconproperties.python_files.getWeightingPotential import get_weighting_potential

    # Planar sensor:

    thickness = 200  # [um]
    width = 40  # [um]

    y, x = np.mgrid[0:thickness:100j, -width * 2:width * 2:100j]
    y_1d = y.T[0]

    # Plot planar electric field in 1 dim
    for V_bias in [50, 150]:
        for n_eff in [1e12, 5e12]:
            E_x, E_y = get_field(x=None,
                                 y=y_1d,
                                 V_bias=V_bias,
                                 n_eff=n_eff,
                                 D=thickness,
                                 is_planar=True)
            plt.plot(y_1d, E_y, linewidth=2.,
                     label='$\mathrm{V_{bias} = %d}$ V, $\mathrm{N_{eff} = %1.2e}$' % (V_bias, n_eff))

    plt.title('Electric field between planar sensor electrodes')
    plt.xlabel('Position [um]')
    plt.ylabel('Electric field [V/um]')
    plt.legend(loc=0)
    plt.grid()
    plt.savefig('Electric_field_planar_1D.pdf', layout='tight')
    plt.show()

    # Plot planar electric field in 2 dim, full fill factor assumed
    E_x, E_y = get_field(x, y, V_bias=20., n_eff=1e12, D=thickness, S=None, is_planar=True)
    plt.quiver(x, y, E_x, E_y, pivot='mid', color='gray', scale=8.)
    plt.savefig('Electric_field_planar.pdf', layout='tight')
    plt.show()

    # 3D sensor:

    radius = 6  # [um]
    distance = 40  # [um]

    y, x = np.mgrid[-distance:distance:200j, -1.618 * distance:1.618 * distance:300j]
    phi_w = get_weighting_potential(x, y, D=distance, S=radius, is_planar=False)

    # Plot weighting potential with colour and contour lines
    levels = np.arange(0, 1, 5)

    plt.contour(x, y, phi_w, 15, colors='black')

    # Plot electrodes
    circle1 = plt.Circle((-distance / 2., 0), radius=radius, color='gray', linewidth=2)
    circle2 = plt.Circle((distance / 2., 0), radius=radius, color='gray', linewidth=2)

    plt.gca().add_patch(circle1)
    plt.gca().add_patch(circle2)

    # Plot electric field in 2D, no neigbouring pixels and space charge assumed
    y, x = np.mgrid[-distance:distance:20j, -1.618 * distance:1.618 * distance:30j]
    E_x, E_y = get_field(x, y, V_bias=20., n_eff=1e12, D=distance, S=radius, is_planar=False)

    plt.quiver(x, y, E_x, E_y, pivot='mid', color='gray', scale=18.)

    plt.title('Electric field (3D sensor)')
    plt.xlabel('Position [um]')
    plt.ylabel('Depth [um]')
    plt.savefig('Electric_field_3D.pdf', layout='tight')
    plt.gca().set_aspect(1.)
    plt.show()
