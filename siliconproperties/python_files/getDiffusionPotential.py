import numpy as np
from scipy import constants


def get_diffusion_potential(n_eff, temperature):
    ''' Diffusion potential [V] as a function of the effective doping
        concentration n_eff [10^12 / cm^3] and the temperature [K].

        Check/citetation of formulars needed!'''

    # [10^12/cm^3] empirical fit at K = 300 range
    N_i = 9.38 * 10 ** 7. * \
        (temperature / 300.)**2 * np.exp(-6884. / temperature)

    # Diffusion potential in thermal equilibrium [V]
    return constants.Boltzmann * temperature / constants.elementary_charge * np.log(n_eff**2. / N_i**2)


if __name__ == '__main__':
    import matplotlib.pylab as plt

    n_eff = np.linspace(0.1, 100., 1000.)

    # Plot diffusion potential
    for temperature in [200, 250, 300, 350]:
        plt.plot(n_eff, get_diffusion_potential(
            n_eff, temperature=temperature), linewidth=2.,
            label='T = %d' % temperature)
    plt.title(
        'Diffusion potential at thermal equilibrium in silicon')
    plt.xlabel('Effective doping concentration [$\mathrm{10^{12} / cm^3}}$]')
    plt.ylabel('Diffusion potential [$\mathrm{V}$]')
    plt.legend(loc=0)
    plt.grid()
    plt.savefig('DiffusionPotential.pdf', layout='tight')
    plt.show()
