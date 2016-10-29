import numpy as np
from scipy import constants

from siliconproperties.python_files.constant import epsilon_s


def get_depletion_voltage(n_eff, distance):
    """ This function returns the full depletion Voltage [V] as a function
        of the effective doping concentration Neff [10^12 /cm^-3] and the
        distance between electrodes [um]. Formular is standard and for example
        used in:
        G. Kramberger et al., Nucl. Inst. And. Meth. A 476 (2002) 645-651
        'Determination of effective trapping times for electrons and holes
        in irradiated silicon'.
    """

    # Constants
    # Relative permitivity of silicon
    epsilon_r = epsilon_s / constants.epsilon_0

    return constants.elementary_charge * n_eff / (constants.epsilon_0 * epsilon_r) * distance ** 2. / 2 * 1e6

if __name__ == '__main__':
    import matplotlib.pylab as plt

    n_eff = np.linspace(0.1, 100., 1000.)

    # Plot depletion voltage
    for distance in [100, 150, 200, 250]:
        plt.plot(n_eff, get_depletion_voltage(
            n_eff, distance=distance), linewidth=2.,
            label='Electrode distance = %d um' % distance)
    plt.title(
        'Full depletion voltage in silicon')
    plt.xlabel('Effective doping concentration [$\mathrm{10^{12} / cm^3}}$]')
    plt.ylabel('Depletion Voltage [$\mathrm{V}}$]')
    plt.legend(loc=0)
    plt.grid()
    plt.savefig('DepletionVoltage.pdf', layout='tight')
    plt.show()
