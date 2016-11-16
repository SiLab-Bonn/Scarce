import matplotlib.pyplot as plt
from fipy.tools import numerix
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.patches import Rectangle
from matplotlib import colors, cm


def get_mesh_plot(fig, mesh, values=None, invert_y_axis=True):
    ax = fig.add_subplot(111)

    vertexIDs = mesh._orderedCellVertexIDs
    vertexCoords = mesh.vertexCoords
    xCoords = numerix.take(vertexCoords[0], vertexIDs)
    yCoords = numerix.take(vertexCoords[1], vertexIDs)

    polys = []

    for x, y in zip(xCoords.swapaxes(0, 1), yCoords.swapaxes(0, 1)):
        if hasattr(x, 'mask'):
            x = x.compressed()
        if hasattr(y, 'mask'):
            y = y.compressed()
        polys.append(zip(x, y))

    collection = PolyCollection(polys)
    collection.set_facecolors('None')
    if invert_y_axis:
        ax.invert_yaxis()
    ax.add_collection(collection)

    if values:
        pass
        #     rgba = cm.cmap(plt.colors.norm(Z))
#         rgba = cm.get_cmap('coolwarm')(colors.norm(Z))
#         collection.set_facecolors(rgba)
#         collection.set_edgecolors(rgba)

    ax.plot()


def plot_mesh(mesh, values=None, invert_y_axis=True):
    fig = plt.figure()
    get_mesh_plot(fig, mesh, values=values, invert_y_axis=invert_y_axis)
    fig.gca().set_aspect('equal')
    plt.show()


def plot_planar_sensor(width,
                       pitch,
                       thickness,
                       n_pixel,
                       V_backplane,
                       V_readout,
                       potential_function=None,
                       field_function=None,
                       mesh=None):
    fig = plt.figure()
    get_planar_sensor_plot(fig,
                           width,
                           pitch,
                           thickness,
                           n_pixel,
                           V_backplane,
                           V_readout,
                           potential_function=potential_function,
                           field_function=field_function,
                           mesh=mesh)
    plt.show()


def get_planar_sensor_plot(fig,
                           width,
                           pitch,
                           thickness,
                           n_pixel,
                           V_backplane,
                           V_readout,
                           potential_function=None,
                           field_function=None,
                           mesh=None):

    ax = fig.add_subplot(111)

    min_x, max_x = -width * float(n_pixel) / 2., width * float(n_pixel) / 2.
    min_y, max_y = 0., thickness

    # Define plot space
    x = np.linspace(min_x, max_x, 1000)
    y = np.linspace(min_y, max_y, 1000)
#     y = np.square(np.linspace(np.sqrt(min_y), np.sqrt(max_y), 1000)) TODO: make this work with streamplot

    # Create x,y plot grid
    xx, yy = np.meshgrid(x, y, sparse=True)

    if potential_function:
        # Plot Potential
        phi = potential_function(xx, yy)
        # BUG in matplotlib: aspect to be set to equal, otherwise contour plot wrong aspect ratio
        # http://stackoverflow.com/questions/28857673/wrong-aspect-ratio-for-contour-plot-with-python-matplotlib
        ax.set_aspect('equal')
        ax.contour(x, y, phi, 10, colors='black')
        cmesh = ax.pcolormesh(x, y, phi, cmap=cm.get_cmap('Blues'), vmin=V_backplane, vmax=V_readout)
        cax = fig.add_axes([ax.get_position().xmax, 0.1, 0.05, ax.get_position().ymax - ax.get_position().ymin])
        fig.colorbar(cmesh, cax=cax, orientation='vertical')

    # Plot E-Field
    if field_function:
        E_x, E_y = field_function(xx, yy)
        ax.streamplot(x, y, E_x, E_y, density=1.0, color='gray', arrowstyle='-')
    elif potential_function:
        E_y, E_x = np.gradient(-phi, np.gradient(y), np.gradient(x))
        ax.streamplot(x, y, E_x, E_y, density=1.0, color='gray', arrowstyle='-')

    if mesh:
        get_mesh_plot(fig, mesh, invert_y_axis=False)

    # Plot backside
    ax.add_patch(Rectangle((min_x, max_y), (max_x - min_x), max_y, color="grey", linewidth=2))

    # Plot pixel(s)
    for pixel in range(n_pixel):
        pixel_position = width * (pixel + 1. / 2.) - width * n_pixel / 2.
        ax.add_patch(Rectangle((pixel_position - pitch / 2, ax.get_ylim()[0]), pitch, 0, color="darkred", linewidth=5))
        ax.plot([pixel_position - pitch / 2, pixel_position - pitch / 2], [min_y, max_y], '--', color='black', linewidth=4)
    ax.plot([pixel_position + pitch / 2, pixel_position + pitch / 2], [min_y, max_y], '--', color='black', linewidth=4)  # Last pixel end

    ax.set_ylim((- 0.02 * (ax.get_ylim()[1] - ax.get_ylim()[0]), 1.02 * max_y))
    ax.set_xlabel('Position x/y [um]', fontsize=22)
    ax.set_ylabel('Position z [um]', fontsize=22)
    ax.invert_yaxis()


def plot_3D_sensor(potential_function, pitch_x, pitch_y, n_pixel, radius, V_readout, V_bias,
                   min_x,
                   max_x,
                   min_y,
                   max_y, nD=2):

    # Interpolate potential to a x,y grid
    xnew = np.linspace(min_x, max_x, 1000)
    ynew = np.linspace(min_y, max_y, 1000)
    xnew_plot, ynew_plot = np.meshgrid(xnew, ynew)
    phi = potential_function(xnew_plot, ynew_plot)

    # Plot Potential
    plt.contour(xnew, ynew, phi, 10, colors='black')
    plt.pcolormesh(xnew_plot, ynew_plot, phi, cmap=cm.get_cmap('Blues'), vmin=V_bias, vmax=V_readout)
    plt.colorbar()

    # Plot E-Field
    dx, dy = np.diff(xnew)[0], np.diff(ynew)[0]
    E_y, E_x = np.gradient(phi, dx, dy)
    E_x, E_y = -E_x, -E_y
    plt.streamplot(xnew_plot, ynew_plot, E_x, E_y, density=1.0, color='gray', arrowstyle='-')

    # Plot readout pillars
    for pillar in range(nD):
        position = pitch_x / nD * (pillar + 1. / 2.) - pitch_x / 2.
        plt.gca().add_patch(plt.Circle((position, 0.), radius, color="darkred", linewidth=0))

    # Plot bias pillars
    positions = ([- pitch_x / 2., -pitch_y / 2.],
                 [0, -pitch_y / 2.],
                 [pitch_x / 2., -pitch_y / 2.],
                 [- pitch_x / 2., pitch_y / 2.],
                 [0, pitch_y / 2.],
                 [pitch_x / 2., pitch_y / 2.])
    for pos_x, pos_y in positions:
        plt.gca().add_patch(plt.Circle((pos_x, pos_y), radius, color="darkblue", linewidth=0))

    plt.xlim((1.05 * min_x, 1.05 * max_x))
    plt.ylim((1.05 * min_y, 1.05 * max_y))
    plt.xlabel('Position x [um]', fontsize=22)
    plt.ylabel('Position y [um]', fontsize=22)
    plt.gca().set_aspect(1.)
#     plt.savefig('3D.svg', dpi=1, layout='tight')
    plt.show()


if __name__ == '__main__':
    from scarce import fields

    thickness = 200  # [um]
    width = 250  # [um]

    def potential_function(x, y):
        return fields.get_weighting_potential_analytic(x, y, D=thickness, S=width, is_planar=True)

    def field_function(x, y):
        return fields.get_weighting_field_analytic(x, y, D=thickness, S=width, is_planar=True)

    from scarce import geometry
    width = 50
    pitch = 50
    n_pixel = 5
    thickness = 200
    resolution = 50
    V_backplane, V_readout = -1, 0

#     mesh=geometry.mesh_planar_sensor(x=width,
#                                      thickness=thickness,
#                                      resolution=resolution)

    plot_planar_sensor(width=width,
                       pitch=width,
                       thickness=thickness,
                       n_pixel=n_pixel,
                       V_backplane=0,
                       V_readout=1,
                       potential_function=potential_function,
                       field_function=field_function,
                       mesh=None)
