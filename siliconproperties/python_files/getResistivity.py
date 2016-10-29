import numpy as np
from scipy import constants

from siliconproperties.python_files.getMobility import get_mobility


def get_resistivity(n_eff, is_n_type=True, temperature=300, e_field=1e3):
    # Calculate the resitivity from:
    # The effective doping concentration n_eff [10^12 / cm^3]
    # the mobility [cm^2/Vs]
    # for n- and p-type silicon. 
    # 
    # The mobility istself is a
    # function of the temperature [K] and the electric field [V/cm].
    # From http://ecee.colorado.edu/~bart/book/mobility.htm

    # TODO: If you take the mobility[E_field] equation seriously, then there is no constant
    # resitivity since the mobility depends also on the electric field. For low E-Fields <= 1000 V/cm
    # the mobility is independent of the E flied and thus the resistivity. Likely this parameter
    # is always given in low field approximation?! Source needed!

    mobility = get_mobility(e_field, temperature, is_electron=is_n_type)

    return 1. / (constants.e * n_eff * mobility * 1e12)

if __name__ == '__main__':
    import matplotlib.pylab as plt

    n_eff = np.logspace(11., 15., 1000.)

    # Plot trapping rate (1 / s)
    plt.plot(n_eff, get_resistivity(n_eff / 1e12, is_n_type=True), label='n-type')
    plt.plot(n_eff, get_resistivity(n_eff / 1e12, is_n_type=False), label='p-type')

    plt.title('Resistivity of silicon (low  e-field approximation)')
    plt.xlabel('Effective doping concentration [$\mathrm{cm^3}}$]')
    plt.ylabel('Resistivity [$\mathrm{\Omega - cm}$]')
    plt.legend(loc=0)
    plt.xscale('log')
    plt.yscale('log')
    plt.grid()
    plt.savefig('Resistivity.pdf', layout='tight')
    plt.show()