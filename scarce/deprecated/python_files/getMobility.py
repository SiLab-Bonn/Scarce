import numpy as np


def get_mobility(e_field, temperature, is_electron):

    # Calculate the mobility [cm^2/Vs] of charge carriers in silicon from
    # the electrical field (E [V/cm]) the temperature (T [K]) and the charge
    # carrier type (isElectron [0/1] otherwise hole).

    # Formular derived from measured data of high purity silicon and the
    # corresponding fit function parameters are used here.
    # From:
    # C. Jacononi et al., Solid state electronics, 1977, vol 20., p. 87
    # 'A review of some charge transport properties of silicon'
    # Note: the doping concentration is irrelevant for < 10^16/cm^3

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


if __name__ == '__main__':
    import matplotlib.pylab as plt

    e_field = np.logspace(3., 5., 1000.)

    # Plot mobility mu
    mu0_e = get_mobility(e_field, temperature=250., is_electron=True)
    mu1_e = get_mobility(e_field, temperature=300., is_electron=True)
    mu0_h = get_mobility(e_field, temperature=250., is_electron=False)
    mu1_h = get_mobility(e_field, temperature=300., is_electron=False)
    # plt.loglog(e_field, mu0_e, linewidth=2., color='blue', linestyle='--', label='Electrons, T = 250K')
    plt.loglog(e_field, mu1_e, linewidth=2., color='blue', linestyle='-', label='Electrons, T = 300K')
    # plt.loglog(e_field, mu0_h, linewidth=2., color='red', linestyle='--', label='Holes, T = 250K')
    plt.loglog(e_field, mu1_h, linewidth=2., color='red', linestyle='-', label='Holes, T = 300K')
    plt.title('Charge carrier mobility in silicon')
    plt.xlabel('Electric field [$\mathrm{V/cm}$]')
    plt.ylabel('Electron/-hole mobility [$\mathrm{cm^2/Vs}$]')
    plt.legend(loc=0)
    plt.grid()
    plt.savefig('Mobility.pdf', layout='tight')
    plt.show()

    # Plot velocity: v = mu * E
    plt.clf()
    v0_e = get_mobility(e_field, temperature=250., is_electron=True) * e_field
    v1_e = get_mobility(e_field, temperature=300., is_electron=True) * e_field
    v0_h = get_mobility(e_field, temperature=250., is_electron=False) * e_field
    v1_h = get_mobility(e_field, temperature=300., is_electron=False) * e_field
    plt.plot(e_field, v0_e, linewidth=2., color='blue', linestyle='-', label='Electrons, T = 250K')
#     plt.plot(e_field, v1_e, linewidth=2., color='blue', linestyle='-', label='Electrons, T = 300K')
    plt.plot(e_field, v0_h, linewidth=2., color='red', linestyle='-', label='Holes, T = 250K')
#     plt.plot(e_field, v1_h, linewidth=2., color='red', linestyle='-', label='Holes, T = 300K')
    plt.title('Charge carrier velocity in silicon')
    plt.xlabel('Electric field [$\mathrm{V/cm}$]')
    plt.ylabel('Electron/-hole velocity [$\mathrm{cm/s}$]')
    plt.plot([5000, 5000], plt.ylim(), '--', color='black', linewidth=2)
    plt.text(6000, 10000000, r'$\frac{\mathrm{100\ V}}{\mathrm{200 \mu m}}$', fontsize=25)
    plt.plot([50000, 50000], plt.ylim(), '--', color='black', linewidth=2)
    plt.text(51000, 10000000, r'$\frac{\mathrm{1000\ V}}{\mathrm{200 \mu m}}$', fontsize=25)
    plt.xlim((1000, 60000))
    plt.legend(loc=0)
    plt.grid()
    plt.gca().set_aspect(1. / plt.gca().get_data_ratio() / 1.618)
    plt.savefig('Velocity.pdf', layout='tight')
    plt.show()
