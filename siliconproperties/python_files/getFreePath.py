import numpy as np

from siliconproperties.python_files.getMobility import get_mobility
from siliconproperties.python_files.getTrapping import get_trapping


def get_free_path(fluence, e_field, temperature, is_electron=True):

    # Calculate the mean free path in cm of the charge carriers from the trapping
    # propability and the velocity. The trapping propability is a
    # function of the fluence and the velocity is a function
    # of the electric field and the temperature. The electric field itself
    # depends on the electrode geometry and the bias voltage.
    velocity = get_mobility(e_field, temperature, is_electron) * e_field
    trapping_time = get_trapping(fluence, is_electron, paper=1)

    return velocity * trapping_time * 10e-9

if __name__ == '__main__':
    import matplotlib.pylab as plt

    fluence = np.logspace(12., 15., 1000.)

    # Plot trapping rate (1 / s)
    s_e, s_e_error = get_free_path(fluence, e_field=1e6, temperature=250, is_electron=True)
    s_h, s_h_error = get_free_path(fluence, e_field=1e6, temperature=250, is_electron=False)
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
    plt.show()
