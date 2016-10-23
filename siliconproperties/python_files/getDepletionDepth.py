import numpy as np
from scipy import constants

from siliconproperties.python_files.getDiffusionPotential import get_diffusion_potential
from siliconproperties.python_files.constant import epsilon_s


def get_depletion_depth(V_bias, n_eff, temperature):
    ''' Depletion depth [um] of silcon with an effective doping concentration n_eff [10^12 /cm^3]
        at a temperature [K] and reverse bias V_bias [V]. 

        Check/citetation of formulars needed!'''

    # Relative permitivity of silicon
    epsilon_r = epsilon_s / constants.epsilon_0

    V_bi = get_diffusion_potential(n_eff=n_eff, temperature=temperature)
    # depletion depth [m]
    return np.sqrt(
        2. * constants.epsilon_0 * epsilon_r * (V_bi + V_bias) / (constants.elementary_charge * n_eff * 10**6))


if __name__ == '__main__':
    import matplotlib.pylab as plt

    V_bias = np.linspace(0, 100., 1000.)

    # Plot diffusion potential
    for n_eff in [1, 10, 100]:
        for temperature in [300]:
            depletion_depth = get_depletion_depth(
                V_bias=V_bias,
                n_eff=n_eff,
                temperature=temperature)
            plt.plot(V_bias,
                     depletion_depth,
                     linewidth=2.,
                     label='$\mathrm{N_{eff}=%d\cdot10^{12}/cm^3}}$, $T=%d$'
                     % (n_eff, temperature))
    plt.title('Depletion depth in silicon')
    plt.xlabel('Bias voltage [$\mathrm{V}}$]')
    plt.ylabel('Depletion depth [$\mathrm{um}$]')
    plt.legend(loc=0)
    plt.grid()
    plt.savefig('DepletionDepth.pdf', layout='tight')
    plt.show()
