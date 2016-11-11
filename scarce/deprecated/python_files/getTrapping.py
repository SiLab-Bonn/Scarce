import numpy as np


def get_trapping(fluence, is_electron, paper=1):

    # Calculate the trapping time tr (e^-(tr) in ns) of charge carriers in
    # silicon as a function of the fluence. There was also a dependence on
    # the temperature measured, that is omitted here!

    if paper == 1:
        # Newer paper where the finding was that there is no difference
        # in electron / hole trapping propability. Electron trapping is similar
        # to Kramberger et. al. but hole trapping much lower.
        # IEEE TRANSACTIONS ON NUCLEAR SCIENCE, VOL. 51, NO. 6, DECEMBER 2004
        # 'Measurement of Trapping Time Constants in
        # Proton-Irradiated Silicon Pad Detectors'

        if is_electron:
            beta = 5.13e-16  # [cm^2/ns]
            beta_error = 0.16e-16  # [cm^2/ns]
        else:
            beta = 5.04e-16  # [cm^2/ns]
            beta_error = 0.18e-16  # [cm^2/ns]

        tr = 1. / (fluence * beta)
        tr_error = 1. / (fluence * beta ** 2) * beta_error

        return tr

    elif paper == 2:
        # Oldest most cited reference, with irradiation to 2 e 14 only
        # Calculate the trapping time tr (e^-(tr) in ns) of charge carriers in
        # silicon at a temperature of T = 263 K as a function of the
        # fluence (in 10^12 Neq/cm^2). This was measured with a laser and planar
        # silicon sensors with a fluence up to 2*10^14 Neq/cm^2. There was a
        # linear behaviour between fluence and the effective trapping
        # propability measured intepended of the silicon type (oxygenated or
        # not and with different doping concentrations) from:
        # G. Kramberger et al., Nucl. Inst. And. Meth. A 476 (2002) 645-651
        # 'Determination of effective trapping times for electrons and holes
        # in irradiated silicon'

        if is_electron:
            beta = 4.2e-16  # [cm^2/ns]
            beta_error = 0.3e-16  # [cm^2/ns]
        else:
            beta = 6.1e-16  # [cm^2/ns]
            beta_error = 0.3e-16  # [cm^2/ns]

        tr = 1. / (fluence * beta)
        tr_error = 1. / (fluence * beta ** 2) * beta_error

        return tr
    else:
        raise RuntimeError('Unknown paper selected!')


if __name__ == '__main__':
    import matplotlib.pylab as plt

    fluence = np.logspace(12., 15., 1000.)

    # Plot trapping rate (1 / s)
    tr_e = get_trapping(fluence, is_electron=True, paper=1)
    tr_h = get_trapping(fluence, is_electron=False, paper=1)
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
    plt.show()