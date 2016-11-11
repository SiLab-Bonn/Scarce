import numpy as np


def get_weighting_potential(x, y, D, S, is_planar=True):
    """ Planar sensor:
        From Nuclear Instruments and Methods in Physics Research A 535 (2004)
        554-557, with correction from wbar = pi*w/2/D to wbar = pi*w/D with:

        x [um] is the offset from the middle of the electrode
        y [um] the position in the sensor
        D [um] the sensor thickness
        S [um] the eletrode width

        3D sensor:
        Weighting potential for two cylinders with:
        S [um] is the radius
        D [um] distance between columns
    """

    # Wheighting potential for one pixel
    if is_planar:
        xbar = np.pi * x / D
        ybar = np.pi * (y - D) / D
        wbar = np.pi * S / D
        return -1. / np.pi * (np.arctan(np.tan(ybar / 2) * np.tanh((xbar + wbar / 2.) / 2.)) -
                              np.arctan(np.tan(ybar / 2) * np.tanh((xbar - wbar / 2.) / 2.)))
    else:
        R = S
        D = D / 2.  # D is the total distance between the columns
        a = np.sqrt(D * D - R * R)
        Phi_w = 1. / (4 * np.arccosh(D / R)) * np.log(((x - a) ** 2 + y ** 2) / ((x + a) ** 2 + y ** 2)) + 0.5

        # Stability
        Phi_w = np.ma.masked_where(np.sqrt((x + D) * (x + D) + y * y) < R, Phi_w)
        Phi_w = np.ma.masked_where(np.sqrt((x - D) * (x - D) + y * y) < R, Phi_w)
        Phi_w = np.ma.masked_where(Phi_w < 0., Phi_w)
        Phi_w = np.ma.masked_where(Phi_w > 1., Phi_w)

        return Phi_w

if __name__ == '__main__':
    import matplotlib.pylab as plt
    from matplotlib import cm

    from siliconproperties.python_files.getWeightingField import (
        get_weighting_field)

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
    plt.savefig('WeightingPotential_planar.pdf', layout='tight')
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
    plt.savefig('WeightingPotential_3D.pdf', layout='tight')
    plt.gca().set_aspect(1.)
    plt.show()
