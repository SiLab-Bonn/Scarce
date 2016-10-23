import numpy as np


def get_weighting_field_planar(x, y, D, S):
    """ From Nuclear Instruments and Methods in Physics Research A 535 (2004)
        554?557, with correction from wbar = pi*w/2/D to wbar = pi*w/D
        with x [um] is the position in the sensor, y [um] the offset from the
        middle of the electrode, D [um] the sensor thickness and S [um] the
        eletrode width. The field is calculated from the drivation of the
        potential in x and y.
    """

    y = D - y  # Electrode at D not at 0
    xbar = np. pi * x / D
    ybar = np.pi * (y - D) / D
    wbar = np.pi * S / D
#     E_x = -1/pi.*( pi.*sin(ybar)./(2.*D.*(cosh(wbar./2+xbar)+cos(ybar))) -  pi.*sin(ybar)./(2.*D.*(cos(ybar)+cosh(wbar./2-xbar))) );    %not easy to find a more simple form
# E_y = -1/pi.*(
# pi.*sinh(wbar./2+xbar)./(2.*D.*(cosh(wbar./2+xbar)+cos(1-pi.*y./D))) +
# pi.*sinh(wbar./2-xbar)./(2.*D.*(cos(1-pi.*y./D)+cosh(wbar./2-xbar))) );
# %not easy to find a more simple form

    E_x = -np.sin(ybar) / (2 * D) * (1. / (np.cosh(xbar - wbar / 2) +
                                           np.cos(ybar)) - 1. / (np.cosh(wbar / 2 + xbar) + np.cos(ybar)))
    E_y = -1. / (2 * D) * (np.sinh(wbar / 2 - xbar) / (np.cosh(wbar / 2 - xbar) +
                                                       np.cos(ybar)) + np.sinh(wbar / 2 + xbar) / (np.cosh(wbar / 2 + xbar) + np.cos(ybar)))

    return E_x, E_y

if __name__ == '__main__':
    import matplotlib.pylab as plt

    thickness = 200  # [um]

    y, x = np.mgrid[0:thickness:100j, -300:300:100j]
    E_x, E_y = get_weighting_field_planar(x, y, D=thickness, S=40)
    
    fig0, ax0 = plt.subplots()
    strm = ax0.streamplot(x, y, E_x, E_y, color=E_x, linewidth=2, cmap=plt.cm.autumn)
    fig0.colorbar(strm.lines)
    
#     fig1, (ax1, ax2) = plt.subplots(ncols=2)
#     ax1.streamplot(x, y, E_x, E_y, density=[0.5, 1])
    
#     lw = 5*speed / speed.max()
#     ax2.streamplot(X, Y, U, V, density=0.6, color='k', linewidth=lw)
    
    plt.show()
# 
#     # Plot diffusion potential
#     for temperature in [200, 250, 300, 350]:
#         plt.plot(n_eff, get_diffusion_potential(
#             n_eff, temperature=temperature), linewidth=2.,
#             label='T = %d' % temperature)
#     plt.title(
#         'Diffusion potential at thermal equilibrium in silicon')
#     plt.xlabel('Effective doping concentration [$\mathrm{10^{12} / cm^3}}$]')
#     plt.ylabel('Diffusion potential [$\mathrm{V}$]')
#     plt.legend(loc=0)
#     plt.grid()
#     plt.savefig('DiffusionPotential.pdf', layout='tight')
#     plt.show()
