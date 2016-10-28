import matplotlib.pyplot as plt
from fipy.tools import numerix
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
    plt.gca().axes.add_collection(collection)

    if values:
        #     rgba = cm.cmap(plt.colors.norm(Z))
        rgba = cm.get_cmap('coolwarm')(colors.norm(Z))
        collection.set_facecolors(rgba)
        collection.set_edgecolors(rgba)

    plt.plot()
    plt.show()
