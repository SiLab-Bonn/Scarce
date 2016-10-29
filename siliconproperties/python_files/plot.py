import matplotlib.pyplot as plt
from fipy.tools import numerix
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib import colors, cm


def plot_mesh(mesh, values=None):

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
    #collection.set_linewidth(0.3 * resolution)
    collection.set_facecolors('white')
    plt.gca().set_aspect(1.)
    plt.gca().axes.add_collection(collection)

    if values:
        #     rgba = cm.cmap(plt.colors.norm(Z))
        rgba = cm.get_cmap('coolwarm')(colors.norm(Z))
        collection.set_facecolors(rgba)
        collection.set_edgecolors(rgba)

    plt.plot()
    plt.show()


def plot_planar_sensor(width,
                            pitch,
                            n_pixel,
                            thickness,potential_fit, V_backplane, V_readout, min_x, max_x, min_y, max_y):
    # Plot potential
    xnew = np.linspace(min_x, max_x, 1000)
    ynew = np.linspace(min_y, max_y, 1000)
    phi = potential_fit(xnew, ynew).T
    plt.contour(xnew, ynew, phi, 15, colors='black')
    plt.pcolormesh(xnew, ynew, phi, cmap=cm.get_cmap('Blues'), vmin=V_backplane, vmax=V_readout)
    plt.colorbar()

    # Plot E-Field
    xnew = np.linspace(min_x, max_x, 10)
    ynew = np.linspace(min_y, max_y, 10)
    z_new = potential_fit(xnew, ynew).T
    xnew, ynew = np.meshgrid(xnew, ynew)
    E_y, E_x = np.gradient(z_new)
    E_x, E_y = -E_x, -E_y
    plt.streamplot(xnew, ynew, E_x, E_y, density=1.5, color='gray', arrowstyle='-')

    plt.show()