''' Script to plot all silicon properties defined in Scarce.
'''

import numpy as np
import matplotlib.pylab as plt

from scarce import silicon


def plot_resistivity():
    V_bias = np.linspace(0, 100., 1000.)
    plt.clf()
    # Plot depletion depth for different bias voltages
    for n_eff in [1, 10, 100]:
        for temperature in [300]:
            depletion_depth = silicon.get_depletion_depth(
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

    # Plot depletion depth as a function of the resistivity
    temperature = 300  # [K]
    e_field = 1e3  # [V/cm]
    V_bias = 100  # [V]
    n_eff = np.logspace(0, 3, 1000.)
    resistivity = silicon.get_resistivity(n_eff,
                                          is_n_type=False,
                                          temperature=temperature,
                                          e_field=e_field)
    depletion_depth = silicon.get_depletion_depth(
        V_bias=V_bias,
        n_eff=n_eff,
        temperature=temperature)
    plt.plot(resistivity,
             depletion_depth,
             linewidth=2.,
             label='$\mathrm{p-type, V_{bias}=%d\ V}$, $T=%d$'
             % (V_bias, temperature))
    plt.plot(resistivity,
             0.3 * np.sqrt(V_bias * resistivity),
             '--',
             linewidth=2.,
             label=r'$0.3\sqrt{V_{bias}\ \mathrm{[V]}\cdot\rho\ \mathrm{[\Omega - cm]}}$')

    resistivity = silicon.get_resistivity(n_eff,
                                          is_n_type=True,
                                          temperature=temperature,
                                          e_field=e_field)
    depletion_depth = silicon.get_depletion_depth(
        V_bias=V_bias,
        n_eff=n_eff,
        temperature=temperature)
    plt.plot(resistivity,
             depletion_depth,
             linewidth=2.,
             label='$\mathrm{n-type, V_{bias}=%d\ V}$, $T=%d$'
             % (V_bias, temperature))

    plt.title('Depletion depth in silicon')
    plt.xlabel('Resistivity [$\mathrm{\Omega - cm}$]')
    plt.ylabel('Depletion depth [$\mathrm{um}$]')
    plt.legend(loc=0)
    plt.grid()
    plt.savefig('DepletionDepthResistivity.pdf', layout='tight')


def plot_depletion_voltage():
    n_eff = np.linspace(0.1, 100., 1000.)
    plt.clf()
    # Plot depletion voltage
    for distance in [100, 150, 200, 250]:
        plt.plot(n_eff, silicon.get_depletion_voltage(
            n_eff, distance=distance), linewidth=2.,
            label='Electrode distance = %d um' % distance)
    plt.title(
        'Full depletion voltage in silicon')
    plt.xlabel('Effective doping concentration [$\mathrm{10^{12} / cm^3}}$]')
    plt.ylabel('Depletion Voltage [$\mathrm{V}}$]')
    plt.legend(loc=0)
    plt.grid()
    plt.savefig('DepletionVoltage.pdf', layout='tight')


def plot_diffusion_potential():
    n_eff = np.linspace(0.1, 100., 1000.)
    plt.clf()
    # Plot diffusion potential
    for temperature in [200, 250, 300, 350]:
        plt.plot(n_eff, silicon.get_diffusion_potential(
            n_eff, temperature=temperature), linewidth=2.,
            label='T = %d' % temperature)
    plt.title(
        'Diffusion potential at thermal equilibrium in silicon')
    plt.xlabel('Effective doping concentration [$\mathrm{10^{12} / cm^3}}$]')
    plt.ylabel('Diffusion potential [$\mathrm{V}$]')
    plt.legend(loc=0)
    plt.grid()
    plt.savefig('DiffusionPotential.pdf', layout='tight')


def plot_eff_acceptor_concentration():
    fluence = np.logspace(0., 4., 1000.)
    plt.clf()
    # Plot diffusion potential
    eff_acc_concentration = silicon.get_eff_acceptor_concentration(
        fluence, n_eff_0=1.8, is_ntype=False, is_oxygenated=False)
    plt.plot(fluence, eff_acc_concentration, '-', linewidth=2.,
             label='p-type')
    eff_acc_concentration = silicon.get_eff_acceptor_concentration(
        fluence, n_eff_0=1.8, is_ntype=False, is_oxygenated=True)
    plt.plot(fluence, eff_acc_concentration, '-.', linewidth=2.,
             label='oxygenated p-type')
    eff_acc_concentration = silicon.get_eff_acceptor_concentration(
        fluence, n_eff_0=1.8, is_ntype=True, is_oxygenated=False)
    plt.plot(fluence, eff_acc_concentration, ':', linewidth=2.,
             label='n-type')
    eff_acc_concentration = silicon.get_eff_acceptor_concentration(
        fluence, n_eff_0=1.8, is_ntype=True, is_oxygenated=True)
    plt.plot(fluence, eff_acc_concentration, '--', linewidth=2.,
             label='oxygenated n-type')

    plt.title('Effective acceptor concentration')
    plt.xlabel('Fluence [$\mathrm{10^{12}\ N_{eq}/cm^2}}$]')
    plt.ylabel('Acceptor concentration $\mathrm{N_{eff}}$ [$\mathrm{10^{12} cm^{-3}}$]')
    plt.legend(loc=0)
    plt.xscale('log')
    plt.yscale('log')
    plt.grid()
    plt.savefig('EffectiveAcceptorConcetration.pdf', layout='tight')


def plot_get_free_path():

    fluence = np.logspace(12., 15., 1000.)
    plt.clf()
    # Plot trapping rate (1 / s)
    s_e = silicon.get_free_path(fluence, e_field=1e6, temperature=250, is_electron=True)
    s_h = silicon.get_free_path(fluence, e_field=1e6, temperature=250, is_electron=False)
    plt.plot(fluence, s_e, linewidth=2., color='blue', linestyle='-', label='Electrons, T = 250')
    plt.plot(fluence, s_h, linewidth=2., color='red', linestyle='-', label='Holes, T = 250')
    plt.title('Charge carrier mean free path in irradiated silicon\nat saturation velocity ($\mathrm{E=10^6\ V/cm}$)')
    plt.xlabel('Fluence [$\mathrm{N_{eq}/cm^2}}$]')
    plt.ylabel('Trapping time [$\mathrm{ns}$]')
    plt.legend(loc=0)
    plt.xscale('log')
    plt.yscale('log')
    plt.grid()
    plt.savefig('MeanFreePath.pdf', layout='tight')


def plot_get_mobility():

    e_field = np.logspace(3., 5., 1000.)
    plt.clf()
    # Plot mobility mu
    mu0_e = silicon.get_mobility(e_field, temperature=250., is_electron=True)
    mu1_e = silicon.get_mobility(e_field, temperature=300., is_electron=True)
    mu0_h = silicon.get_mobility(e_field, temperature=250., is_electron=False)
    mu1_h = silicon.get_mobility(e_field, temperature=300., is_electron=False)
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

    # Plot velocity: v = mu * E
    plt.clf()
    v0_e = silicon.get_mobility(e_field, temperature=250., is_electron=True) * e_field
    v1_e = silicon.get_mobility(e_field, temperature=300., is_electron=True) * e_field
    v0_h = silicon.get_mobility(e_field, temperature=250., is_electron=False) * e_field
    v1_h = silicon.get_mobility(e_field, temperature=300., is_electron=False) * e_field
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


def plot_get_resistivity():
    n_eff = np.logspace(11., 15., 1000.)
    plt.clf()
    # Plot trapping rate (1 / s)
    plt.plot(n_eff, silicon.get_resistivity(n_eff / 1e12, is_n_type=True), label='n-type')
    plt.plot(n_eff, silicon.get_resistivity(n_eff / 1e12, is_n_type=False), label='p-type')

    plt.title('Resistivity of silicon (low  e-field approximation)')
    plt.xlabel('Effective doping concentration [$\mathrm{cm^3}}$]')
    plt.ylabel('Resistivity [$\mathrm{\Omega - cm}$]')
    plt.legend(loc=0)
    plt.xscale('log')
    plt.yscale('log')
    plt.grid()
    plt.savefig('Resistivity.pdf', layout='tight')


def plot_trapping():

    fluence = np.logspace(12., 15., 1000.)

    # Plot trapping rate (1 / s)
    tr_e = silicon.get_trapping(fluence, is_electron=True, paper=1)
    tr_h = silicon.get_trapping(fluence, is_electron=False, paper=1)

    plt.clf()
    plt.plot(fluence, tr_e, linewidth=2., color='blue', linestyle='-', label='Electrons')
    plt.plot(fluence, tr_h, linewidth=2., color='red', linestyle='-', label='Holes')
    plt.title('Charge carrier trapping time in irradiated silicon')
    plt.xlabel('Fluence [\mathrm{$N_{eq}/cm^2}}$]')
    plt.ylabel('Trapping time [$\mathrm{ns}$]')
    plt.legend(loc=0)
    plt.xscale('log')
    plt.yscale('log')
    plt.grid()
    plt.savefig('TrappingTime.pdf', layout='tight')
    
def create_plots():
    plot_resistivity()
    plot_depletion_voltage()
    plot_diffusion_potential()
    plot_eff_acceptor_concentration()
    plot_get_free_path()
    plot_get_mobility()
    plot_get_resistivity()
    plot_trapping()

if __name__ == '__main__':
    create_plots()