import numpy as np

from siliconproperties.python_files.getDepletionVoltage import get_depletion_voltage


def get_potential_planar(x, V_bias, n_eff, D):
    """ Calculates the potential [V] in a planar sensor as a
        function of the position x between the electrodes [um],
        the bias voltage V_bias [V], the effective doping
        concentration n_eff [cm^-3] and the sensor Width D [um].
        The analytical function from the detector book p. 93 is used.
    """

    V_dep = get_depletion_voltage(n_eff, D)  # Depletion voltage

    a = (V_bias - V_dep) / D
    b = -2. * V_dep / (D**2)
    return (a - b / 2 * x) * x

if __name__ == '__main__':
    import matplotlib.pylab as plt

    thickness = 200

    # Plot planar potential
    x = np.linspace(0., thickness, 1000.)

    for V_bias in [50, 150]:
        for n_eff in [1e12, 5e12]:
            phi = get_potential_planar(
                x, V_bias=V_bias, n_eff=n_eff, D=thickness)
            plt.plot(x, phi, linewidth=2.,
                     label='$\mathrm{V_{bias} = %d}$ V, $\mathrm{N_{eff} = %1.2e}$' % (V_bias, n_eff))

    plt.title('Potential between planar sensor electrodes')
    plt.xlabel('Position [um]')
    plt.ylabel('Potential [V]')
    plt.legend(loc=0)
    plt.grid()
    plt.savefig('PlanarPotential.pdf', layout='tight')
    plt.show()
