import numpy as np
from scipy import constants


def getdedx(betagamma, thickness, density_correction=True, restricted=True):
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
    re = 2.817940325 * 10**-13 	# electron radius [cm]
    Z = 14                      # atomic number [#nucleons]
    A = 28.0855                 # atomic mass of silicon absorber [amu = g/mol]
    Na = 6.0221415 * 10**23        # avogrado number
    me = 0.510998918            # mass of an electron [MeV/c**2]
    I = 173 * 10**(-6)             # mean excitation energy [MeV]
    K = 4 * np.pi * Na * re**2 * me         # constant(?) [MeV]
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
    if density_correction and restricted:
        dedx = (-K * z ** 2 * Z) / (A * beta2) * (0.5 * np.log((2 * me * betagamma2 * T_cut) / I ** 2) - beta2 /
                                                  2 * (1 + T_cut / Tmax) - delta / 2)  # dEdx with density effect and restriction on measured dEdx in thin absorbers
        return dedx

#     # dedx with MP Eloss correction
#     kappa = K / 2 * Z / A * mat_length / beta2
#     dedx_density_MPV = kappa * \
#         (np.log(2 * me * betagamma2 / I) +
#          np.log(kappa / I) + j - beta2 - delta) / (mat_length)

#     # MP charge deposit
#     # Mean electron hole pairs measured, FIXME (temp dependent)
#     electron_hole_pair = 3.65  # energy neede to create electron hole pair [eV]
#     dedx_Q = dedx_density_MPV / electron_hole_pair * \
#         1e3 * mat_length  # factor 1e3 for scaling to ke

    return dedx

if __name__ == '__main__':
    import matplotlib.pylab as plt

    # Plot dEdX
    betagamma = np.logspace(-1, 3.5, 1000)
    dedx = getdedx(betagamma=betagamma, thickness=200, density_correction=False, restricted=False)
    dedx_den = getdedx(betagamma=betagamma, thickness=200, density_correction=True, restricted=False)
    dedx_den_res = getdedx(betagamma=betagamma, thickness=200, density_correction=True, restricted=True)

    plt.xscale('log')
    plt.plot(betagamma, -dedx, label='Bethe Bloch', linewidth=2)
    plt.plot(betagamma, -dedx_den, label='+ density correction', linewidth=2)
    plt.plot(betagamma, -dedx_den_res, label='+ restricted ELoss', linewidth=2)
    plt.xlabel(r'$\mathrm{\beta\gamma} = p/(Mc)}$')
    plt.ylabel(r'Average stopping power [$\mathrm{MeV / g\ cm^2}$]')
    plt.title('Average stopping power (energy loss per material distance)')
    plt.grid()
    plt.legend(loc=0)
    plt.show()
#     h = semilogx(betagamma,-dedx, 'LineWidth', 2)
#     hold on;
#     semilogx(betagamma,-dedx_density, 'LineWidth', 2, 'color', 'red');
#     semilogx(betagamma,-dedx_density_rest, 'LineWidth', 2, 'color', 'green');
#     semilogx(betagamma,dedx_density_MPV, 'LineWidth', 2, 'color', 'black');
#     #semilogx(betagamma,delta, 'LineWidth', 2, 'color', 'black');
#     hold off;
#     set(gca, 'YMinorTick','on');
#     title('Average stopping power (energy loss per material distance) ', 'FontWeight','bold','FontSize', 10);
#     xlabel('\beta\gamma = p/(Mc)', 'FontWeight','bold');
#     ylabel('average stopping power [MeV g**{-1} cm**2]', 'FontWeight','bold');
#     xlim([0.3 2000]);
#     ylim([0 5]);
#     set(gcf, 'Color', [1 1 1]);
#     hleg1 = legend('Bethe Bloch','+ density correction','+ restricted Eloss','+ MP Eloss');
#     grid on
# #     xlabh = get(gca,'XLabel');
# #     x2 = axes('Position',get(gca,'Position')-[0.0 0.05 0 0]);
# #     set(x2,'ytick',[1 2 3 4 5], 'yticklabel',{})
#
#     #axes('Position',[0 0 1 1]);
#     #export_fig('dEdx.pdf');
#     figure;
#     semilogx(betagamma,dedx_Q, 'LineWidth', 2);
#     min(dedx_Q)
#     set(gca, 'YMinorTick','on');
#     title_str = sprintf('Most propable charge deposited in #d um silicon layer', thickness);
#     title(title_str, 'FontWeight','bold','FontSize', 10);
#     xlabel('\beta\gamma = p/(Mc)', 'FontWeight','bold');
#     ylabel('most propable charge deposited [ke]', 'FontWeight','bold');
#     set(gca, 'YMinorTick','on');
#     xlim([0.55 2000]);
#     ylim([0 60]);
#     grid on
#     set(gcf, 'Color', [1 1 1]);
#     export_fig('ChargeDeposited.pdf');
#     min(dedx_Q)
