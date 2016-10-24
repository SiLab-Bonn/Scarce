import numpy as np


def get_weighting_field(x, y, D, S, is_planar=True):
    """ From Nuclear Instruments and Methods in Physics Research A 535 (2004)
        554-557, with correction from wbar = pi*w/2/D to wbar = pi*w/D
        with x [um] is the position in the sensor, y [um] the offset from the
        middle of the electrode, D [um] the sensor thickness and S [um] the
        eletrode width. The field is calculated from the drivation of the
        potential in x and y.
    """

    if is_planar:
        # y = D - y  # Electrode at D not at 0
        xbar = np. pi * x / D
        ybar = np.pi * (y - D) / D
        wbar = np.pi * S / D

        # Not easy to find a more simple form
        E_x = -np.sin(ybar) / (2. * D) * (1. / (np.cosh(xbar - wbar / 2.) +
                                                np.cos(ybar)) - 1. / (np.cosh(wbar / 2 + xbar) + np.cos(ybar)))
        E_y = -1. / (2 * D) * (np.sinh(wbar / 2 - xbar) / (np.cosh(wbar / 2 - xbar) +
                                                           np.cos(ybar)) + np.sinh(wbar / 2 + xbar) / (np.cosh(wbar / 2 + xbar) + np.cos(ybar)))

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

        return E_x, E_y

if __name__ == '__main__':
    import matplotlib.pylab as plt
    from matplotlib import cm

    from siliconproperties.python_files.getWeightingPotential import (
        get_weighting_potential)

    thickness = 200  # [um]
    width = 40  # [um]

    y, x = np.mgrid[0:thickness:100j, -width * 2:width * 2:100j]
    phi_w = get_weighting_potential(x, y, D=thickness, S=width, is_planar=True)

    # Plot weighting potential with colour and contour lines
    levels = np.arange(0, 1, 5)
    plt.imshow(phi_w, extent=[x.min(), x.max(), y.max(), y.min()], origin='upper', cmap=cm.get_cmap('Blues'))
    plt.contour(x, y, phi_w, 15, colors='black')

    # Plot electrode
    plt.plot([-width / 2., width / 2.], [0, 0], linewidth=5, color='red')

    # Plot weighting field directions
    y, x = np.mgrid[0:thickness:20j, -width * 2:width * 2:20j]
    E_x, E_y = get_weighting_field(x, y, D=thickness, S=width, is_planar=True)

    plt.quiver(x, y, E_x / np.sqrt(E_x ** 2 + E_y ** 2), E_y / np.sqrt(E_x ** 2 + E_y ** 2), pivot='mid', color='gray', scale=30.)

    plt.title('Weighting potential and weighting field direction (planar sensor)')
    plt.xlabel('Position [um]')
    plt.ylabel('Depth [um]')
    plt.gca().set_aspect(1. / plt.gca().get_data_ratio() / 1.618)
    plt.savefig('WeightingField_planar.pdf', layout='tight')
    plt.show()

    radius = 6  # [um]
    distance = 40  # [um]

    y, x = np.mgrid[-distance:distance:200j, -1.618 * distance:1.618 * distance:300j]
    phi_w = get_weighting_potential(x, y, D=distance, S=radius, is_planar=False)

    # Plot weighting potential with colour and contour lines
    levels = np.arange(0, 1, 5)

#     plt.imshow(phi_w, extent=[x.min(), x.max(), y.max(), y.min()], origin='upper', cmap=cm.get_cmap('Blues'))
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
    plt.xlabel('Position [um]')
    plt.ylabel('Depth [um]')
    plt.savefig('WeightingField_3D.pdf', layout='tight')
    plt.gca().set_aspect(1.)
    plt.show()
