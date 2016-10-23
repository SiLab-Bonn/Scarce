import numpy as np


def get_eff_acceptor_concentration(fluence, n_eff_0, is_ntype, is_oxygenated):
    """ Calculate the effective acceptor concentration [10^12 cm^-3] of irradiated n and
        p-type silicon with and without oxygen enriched as a function of the
        fluence [10^12 cm^-2]. The data can be desribed by different line fits.
        The parameters were extracted from a plot taken for n-type silicon from:
        CERN/LHCC 2000-009 LEB Status Report/RD48 31 Dec. 1999
        and for p-type silicon from:
        RD50 Status Report 2006 CERN-LHCC-2007-005, p. 4-6

        Due to the difference in the data for different technologies
        a rather large error on the propotionality factor of 10% is assumed.
    """

    # Constants
    if is_ntype:  # n-type silicon
        # [cm^-1], slope before type inversion
        beta_n = -5.5e-2
        # [cm^-1], slope after type inversion, standart silicon FZ
        beta_n_fz = 2.83e-2
        # [cm^-1], slope after type inversion, oxygenated silicon FZ
        beta_n_ofz = 0.94e-2
        Phi_inv = -n_eff_0 / beta_n  # Inversion point N_eq [10^14 cm^-2]

        # After type inversion
        if is_oxygenated:  # Oxygenated silicon
            n_eff = beta_n_ofz * (fluence - Phi_inv)
        else:  # Standart silicon
            n_eff = beta_n_fz * (fluence - Phi_inv)

        # Before type inversion
        n_eff[fluence < Phi_inv] = n_eff_0 + \
            beta_n * fluence[fluence < Phi_inv]

    else:  # p-type silicons
        beta_n_fz = 2.3e-2  # [cm^-1]
        beta_n_ofz = 2.7e-2  # [cm^-1]
        if is_oxygenated:  # Oxygenated silicon
            n_eff = n_eff_0 + beta_n_ofz * fluence
        else:  # Standart silicon
            n_eff = n_eff_0 + beta_n_fz * fluence

    n_eff_error = n_eff * 0.1

    return n_eff, n_eff_error

if __name__ == '__main__':
    import matplotlib.pylab as plt

    fluence = np.logspace(12., 15., 1000.)

    # Plot diffusion potential
    eff_acc_concentration = get_eff_acceptor_concentration(
        fluence, n_eff_0=1.8e12, is_ntype=False, is_oxygenated=False)[0]
    plt.plot(fluence, eff_acc_concentration, linewidth=2.,
             label='p-type')
    eff_acc_concentration = get_eff_acceptor_concentration(
        fluence, n_eff_0=1.8e12, is_ntype=False, is_oxygenated=True)[0]
    plt.plot(fluence, eff_acc_concentration, linewidth=2.,
             label='oxygenated p-type')
    eff_acc_concentration = get_eff_acceptor_concentration(
        fluence, n_eff_0=1.8e12, is_ntype=True, is_oxygenated=False)[0]
    plt.plot(fluence, eff_acc_concentration, linewidth=2.,
             label='n-type')
    eff_acc_concentration = get_eff_acceptor_concentration(
        fluence, n_eff_0=1.8e12, is_ntype=True, is_oxygenated=True)[0]
    plt.plot(fluence, eff_acc_concentration, linewidth=2.,
             label='oxygenated n-type')

    plt.title('Effctive acceptor concentration')
    plt.xlabel('Fluence [$\mathrm{N_{eq}/cm^2}}$]')
    plt.ylabel('Effective acceptor concentration [$\mathrm{N_{eff}\ 10^{12} cm^{-3}}$]')
    plt.legend(loc=0)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.grid()
    plt.savefig('EffectiveAcceptorConcetration.pdf', layout='tight')
    plt.show()
