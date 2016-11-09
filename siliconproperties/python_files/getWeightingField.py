import numpy as np


def get_weighting_field(x, y, D, S, is_planar=True):
    """ From Nuclear Instruments and Methods in Physics Research A 535 (2004)
        554-557, with correction from wbar = pi*w/2/D to wbar = pi*w/D
        with x [um] is the position in the sensor [0:thickness], y [um] the offset from the
        middle of the electrode, D [um] the sensor thickness and S [um] the
        eletrode width. The field is calculated from the drivation of the
        potential in x and y.
    """

    if is_planar:

#         y = D - y  # Electrode at D not at 0
        xbar = np.pi * x / D
        ybar = np.pi * (y) / D
        wbar = np.pi * S / D

#         # Likely no simpler form possible
#         E_x = np.sin(ybar) / (2. * D) * (1. / (np.cosh(xbar - wbar / 2.) + np.cos(ybar)) 
#                                         - 1. / (np.cosh(xbar + wbar / 2.) + np.cos(ybar)))
#         
#         E_y = 1. / (2 * D) * (-np.sinh(- xbar + wbar / 2.) / (np.cosh(- xbar + wbar / 2.) + np.cos(ybar)) 
#                               -np.sinh(  xbar + wbar / 2.) / (np.cosh(  xbar + wbar / 2.) + np.cos(ybar)))
         
#         return E_x, E_y
        
        # Not easy to find a more simple form
        
        
        
        denom = (np.cosh(1./2. * (wbar- 2. * xbar)) + np.cos(ybar)) * (np.cosh(1./2. * (wbar + 2. * xbar)) + np.cos(ybar))
        
#         print np.sin(ybar) * np.sinh(wbar/2.) * np.sinh(xbar) / denom- E_x
#         print np.allclose(- np.sinh(wbar/2.) * (np.cos(wbar/2.) + np.cos(ybar)*np.cosh(xbar)) / denom, E_y)

        
        E_x = - np.sin(ybar) * np.sinh(wbar/2.) * np.sinh(xbar) / denom
        
        E_y = np.sinh(wbar/2.) * (np.cosh(wbar/2.) + np.cos(ybar)*np.cosh(xbar)) / denom
       
        return E_x, E_y
        

    else:
        # 3D sensor:
        # From the analytical derivation of the get_weighting_potential function
        # Weighting potential for two cylinders with:
        # S [um] is the radius
        # D [um] distance between columns

        R = S
        D = D / 2.
        a = np.sqrt(D * D - R * R)

        E_x = a / (np.arccosh(D / R)) * (a ** 2 - x ** 2 + y ** 2) / (((a - x) ** 2 + y ** 2) * ((a + x) ** 2 + y ** 2))
        E_y = -2 * a / (np.arccosh(D / R)) * (x * y) / (((a - x) ** 2 + y ** 2) * ((a + x) ** 2 + y ** 2))

        E_x = np.ma.masked_where(np.sqrt((x + D) * (x + D) + y * y) < R, E_x)
        E_x = np.ma.masked_where(np.sqrt((x - D) * (x - D) + y * y) < R, E_x)
        E_y = np.ma.masked_where(np.sqrt((x + D) * (x + D) + y * y) < R, E_y)
        E_y = np.ma.masked_where(np.sqrt((x - D) * (x - D) + y * y) < R, E_y)

        return -E_x, -E_y

if __name__ == '__main__':
    import matplotlib.pylab as plt
    from matplotlib import cm

    from siliconproperties.python_files.getWeightingPotential import (
        get_weighting_potential)

    thickness = 1.  # [um]
    width = 1.  # [um]

    x = np.linspace(-width * 2, width * 2, 200)
    y = np.linspace(0, thickness, 200)
    xx, yy = np.meshgrid(x, y, sparse=True)
#     y_1d = y.T[0]

#     # Plot planar weighting field in 1 dim
#     for pitch, line_style in zip([50, 250, 500], ['-', '-.', '--']):
#         phi_w = get_weighting_potential(np.zeros_like(y_1d), y_1d, D=thickness, S=pitch, is_planar=True)
#         E_x, E_y = get_weighting_field(x=np.zeros_like(y_1d),
#                                        y=y_1d,
#                                        D=thickness,
#                                        S=pitch,
#                                        is_planar=True)
#         plt.plot(y_1d, phi_w, linestyle=line_style, linewidth=2.,
#                  label='$\mathrm{\Phi_{w},\ %d\ \mu m\ pitch}$' % pitch, color='red')
# 
#     plt.ylabel('Weighting potential [$\mathrm{V}$]')
#     plt.grid()
#     plt.xlabel('Depth [$\mathrm{\mu m}$]')
# 
#     h1, l1 = plt.gca().get_legend_handles_labels()
#     ax2 = plt.twinx(plt.gca())
# 
#     # Plot planar weighting potential in 1 dim
#     for pitch, line_style in zip([50, 250, 500], ['-', '-.', '--']):
#         phi_w = get_weighting_potential(np.zeros_like(y_1d), y_1d, D=thickness, S=pitch, is_planar=True)
#         E_x, E_y = get_weighting_field(x=np.zeros_like(y_1d),
#                                        y=y_1d,
#                                        D=thickness,
#                                        S=pitch,
#                                        is_planar=True)
#         ax2.plot(y_1d, E_y, linestyle=line_style, linewidth=2.,
#                  label='$\mathrm{|E|_w,\ %d\ \mu m\ pitch}$' % pitch, color='blue')
# 
#     h2, l2 = ax2.get_legend_handles_labels()

#     plt.title('Weighting potential and field at planar pixel center, $\mathrm{d = 300\ \mu m}$')
#     plt.legend(h1 + h2, l1 + l2, loc=0)
#     plt.ylim((0, 0.05))
#     plt.grid()
#     plt.ylabel('Weighting field [$\mathrm{V/\mu m}$]')
#     plt.xlabel('Depth [$\mathrm{\mu m}$]')
#     plt.savefig('WeightingField_planar_1D.pdf', layout='tight')
#     plt.show()


#     print get_weighting_field(3., 0.1, D=1., S=1., is_planar=True)
#     raise

#     phi_w = get_weighting_potential(xx, yy, D=thickness, S=width, is_planar=True)
#  
#     # Plot weighting potential with colour and contour lines
#     # BUG in matplotlib: aspect to be set first, otherwise contour plot wrong
#     plt.gca().set_aspect('equal')
#     plt.contour(x, y, phi_w, 10, colors='black')
#     plt.pcolormesh(x, y, phi_w, cmap=cm.get_cmap('Blues'), vmin=phi_w.min(), vmax=phi_w.max())
#     plt.colorbar()
#     # Plot electrode
# #     plt.plot([-width / 2., width / 2.], [0, 0], linewidth=5, color='red')
# 
#     # Plot weighting field directions
#     E_x, E_y = get_weighting_field(xx, yy, D=thickness, S=width, is_planar=True)
#     E_x_2, E_y_2 = np.gradient(phi_w, np.diff(x)[0], np.diff(y)[0])
# 
#     plt.streamplot(x, y, E_x_2, E_y_2, density=1.0, color='black', arrowstyle='->')
#     plt.streamplot(x, y, E_x, E_y, density=1.0, color='gray', arrowstyle='->')
#     print plt.gca().get_data_ratio()
# 
#     plt.title('Weighting potential and weighting field direction (planar sensor)')
#     plt.xlabel('Position [um]')
#     plt.ylabel('Depth [um]')
# #     plt.gca().invert_yaxis()
#     plt.contour(x, y, phi_w, 10, colors='white')
# #     plt.savefig('WeightingField_planar.pdf', layout='tight')
#     plt.show()

    plt.clf()
    y = np.linspace(0., 1., 100)
    x = np.zeros_like(y)

    phi_w = get_weighting_potential(x, y, D=1, S=1)
    
    plt.plot(y, phi_w)
    
    plt.show()

    raise

    radius = 6  # [um]
    distance = 70  # [um]

    y, x = np.mgrid[-distance:distance:200j, -1.618 * distance:1.618 * distance:300j]
    x_1d = x[0]

    # Plot planar weighting field in 1 dim
    phi_w = get_weighting_potential(x_1d, np.zeros_like(x_1d), D=distance, S=radius, is_planar=False)
    plt.plot(x_1d, phi_w, linestyle='-', linewidth=2.,
             label='$\mathrm{\Phi_{w}}$', color='red')

    plt.ylabel('Weighting potential [$\mathrm{V}$]')
    plt.grid()
    plt.ylim((-0.2, 1.2))

    # Plot electrode
    plt.gca().add_patch(plt.Rectangle(((-distance - radius) / 2., plt.ylim()[0]), radius, plt.ylim()[1] - plt.ylim()[0], color="grey"))
    plt.gca().add_patch(plt.Rectangle(((distance - radius) / 2., plt.ylim()[0]), radius, plt.ylim()[1] - plt.ylim()[0], color="grey"))
    plt.xlabel('Depth [$\mathrm{\mu m}$]')

    h1, l1 = plt.gca().get_legend_handles_labels()
    ax2 = plt.twinx(plt.gca())

    # Plot planar weighting potential in 1 dim
    E_x, E_y = get_weighting_field(x=x_1d,
                                   y=np.zeros_like(x_1d),
                                   D=distance,
                                   S=radius,
                                   is_planar=False)
    ax2.plot(x_1d, E_x, linestyle='-.', linewidth=2.,
             label='$\mathrm{E_{x, w}}$', color='blue')

    h2, l2 = ax2.get_legend_handles_labels()

    plt.title('Weighting potential and field at 3D pixel center line, $\mathrm{d = %d\ \mu m}$' % distance)
    plt.grid()
    
    plt.legend(h1 + h2, l1 + l2, loc=0)
    plt.ylabel('Weighting field [$\mathrm{V/\mu m}$]')
    plt.xlabel('Position [$\mathrm{\mu m}$]')
    plt.savefig('WeightingField_3D_1D.pdf', layout='tight')
    plt.show()

    phi_w = get_weighting_potential(x, y, D=distance, S=radius, is_planar=False)

    # Plot weighting potential with colour and contour lines
    levels = np.arange(0, 1, 5)

    plt.contour(x, y, phi_w, 15, colors='black')
    plt.pcolor(x, y, phi_w, cmap='Blues', vmin=np.min(phi_w), vmax=np.max(phi_w))

    # Plot electrodes
    circle1 = plt.Circle((-distance / 2., 0), radius=radius, color='gray', linewidth=2)
    circle2 = plt.Circle((distance / 2., 0), radius=radius, color='gray', linewidth=2)

    plt.gca().add_patch(circle1)
    plt.gca().add_patch(circle2)

    # Plot weighting field directions
    y, x = np.mgrid[-distance:distance:20j, -1.618 * distance:1.618 * distance:30j]
    E_x, E_y = get_weighting_field(x, y, D=distance, S=radius, is_planar=False)

    plt.quiver(x, y, E_x / np.sqrt(E_x ** 2 + E_y ** 2), E_y / np.sqrt(E_x ** 2 + E_y ** 2), pivot='mid', color='gray', scale=30.)

    plt.title('Weighting potential and weighting field direction (3D sensor)')
    plt.xlabel('Position X [um]')
    plt.xlabel('Position Y [um]')
    plt.savefig('WeightingField_3D.pdf', layout='tight')
    plt.gca().set_aspect(1.)
    plt.show()
