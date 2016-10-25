import numpy as np
from scipy import interpolate

from siliconproperties.python_files import constant


def get_attenuation(energy, process='total'):
    """ Interpolates the data from http://physics.nist.gov/cgi-bin/Xcom
    to return the attenuation length for the given photon energy [MeV] 
    in silicon for:
        - compton effect
        - photo effect
        - pair production
        - total
    """

    data = np.loadtxt(fname='AttenuationData.txt',
                      dtype=np.float,
                      skiprows=11).T

    x = data[0]

    if 'compton' in process.lower():
        y = data[2]
    elif 'photo' in process.lower():
        y = data[3]
    elif 'pair' in process.lower():
        y = data[4] + data[5]
    elif 'total' in process.lower():
        y = data[6]
    else:
        raise AttributeError('Unknown process %s. Try compton, photo, pair, total.', process)

    y *= constant.density  # cm2 / g -> 1 / cm

    interpol = interpolate.interp1d(x, y, bounds_error=False, fill_value=np.NaN)  # , kind='quadratic')
    y_interpolated = interpol(energy)

    return y_interpolated

if __name__ == '__main__':
    import matplotlib.pylab as plt

    plt.ylim((1e-3, 1e5))

    # Add some common X-Ray ources
    sources_names = point_label = [r'$Fe\ K_{\alpha}$', r'$Mo\ K_{\alpha}$', r'$Ba\ K_{\alpha}$', r'$Am_{2,1}$']
    sources_energies = np.array([5.89875, 17.426, 32.004, 59.541])  # keV
    for index in range(len(sources_energies)):
        source_energy = sources_energies[index] / 1000.
        plt.plot([source_energy, source_energy], plt.ylim(), color='lightgrey', linewidth=2)
        plt.gca().annotate('{}'.format(sources_names[index]), xy=(source_energy, get_attenuation(source_energy) * 15.), xytext=(20, 12), ha='right', textcoords='offset points', fontsize=14)

    # Plot attenuation
    energy = np.logspace(-3, 3, 2000)
    processes = ['Compton effect', 'Photo effect', 'Pair production', 'Total']
    marker_styles = ['--', '-.', ':', '-']
    colors = ['gray'] * 3 + ['black']
    for index, process in enumerate(processes):
        plt.plot(energy, get_attenuation(energy, process=process), marker_styles[index], linewidth=2, label=process, color=colors[index])

    plt.xscale('log')
    plt.yscale('log')
    plt.grid()
    plt.title('Attenuation for photons in silicon (Z=14)')
    plt.xlabel('Photon Energy [MeV]')
    plt.ylabel('Attenutation [1/cm]')
    plt.legend(loc=0)
    plt.savefig('Attenuation.pdf', layout='tight')
    plt.show()
