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
        pass
        #     rgba = cm.cmap(plt.colors.norm(Z))
#         rgba = cm.get_cmap('coolwarm')(colors.norm(Z))
#         collection.set_facecolors(rgba)
#         collection.set_edgecolors(rgba)

    plt.plot()
    plt.show()


def plot_planar_sensor(potential_function,
                       width,
                       pitch,
                       n_pixel,
                       thickness, 
                       V_backplane, 
                       V_readout, 
                       min_x, max_x, 
                       min_y, max_y,
                       field_function=None):

    # Plot potential
    xnew = np.linspace(min_x, max_x, 1000)
    ynew = np.linspace(min_y, max_y, 1000)
    phi = potential_function(xnew, ynew).T
    ynew_plot = np.linspace(min_y, max_y, 1000)
    plt.contour(xnew, ynew, phi, 15, colors='black')
    plt.pcolormesh(xnew, ynew_plot, phi, cmap=cm.get_cmap('Blues'), vmin=V_backplane, vmax=V_readout)
    plt.colorbar()

    # Plot E-Field
    xnew = np.linspace(min_x, max_x, 10)
    ynew = np.linspace(min_y, max_y, 10)
    z_new = potential_function(xnew, ynew).T
    xnew_plot, ynew_plot = np.meshgrid(xnew, ynew)
    E_y, E_x = np.gradient(z_new)
    E_x, E_y = -E_x, -E_y
    plt.streamplot(xnew_plot, ynew_plot, E_x, E_y, density=1.5, color='gray', arrowstyle='-')
    
    # Plot pixel background
    plt.gca().add_patch(plt.Rectangle((min_x, plt.ylim()[1]), (max_x - min_x), 0.05*(plt.ylim()[1] - plt.ylim()[0]), color="grey", linewidth=0))
    
    # Plot pixel(s)
    for pixel in range(n_pixel):
        pixel_position = width * (pixel + 1. / 2.) - width * n_pixel / 2.
        plt.gca().add_patch(plt.Rectangle((pixel_position - pitch / 2, plt.ylim()[1]), pitch, 0.05*(plt.ylim()[1] - plt.ylim()[0]), color="darkblue", linewidth=0))
#         plt.gca().add_patch(plt.Rectangle(((distance - radius) / 2., plt.ylim()[0]), radius, plt.ylim()[1] - plt.ylim()[0], color="grey"))
    
    # Plot backside
    plt.gca().add_patch(plt.Rectangle((min_x, plt.ylim()[0]), (max_x - min_x), - 0.05*(plt.ylim()[1] - plt.ylim()[0]), color="darkred", linewidth=0))
    
    plt.ylim((1.05 * min_y, 1.05 * max_y))
    plt.xlabel('Position [um]')
    plt.ylabel('Depth [um]')
    plt.show()
    

def plot_3D_sensor():
    pass

if __name__ == '__main__':
    from siliconproperties.python_files.getWeightingPotential import (
        get_weighting_potential)
    from siliconproperties.python_files.getWeightingField import (
        get_weighting_field)
    
    thickness = 200  # [um]
    width = 40  # [um]

    y, x = np.mgrid[0:thickness:100j, -width * 2:width * 2:100j]
    
    def potential_function(x, y):
        return get_weighting_potential(x, y, D=thickness, S=width, is_planar=True)
    
    plot_planar_sensor(potential_function=potential_function, 
                       width=width, 
                       pitch=width, 
                       n_pixel=1, 
                       thickness=thickness, 
                       V_backplane=0, 
                       V_readout=1, 
                       min_x=0, 
                       max_x=thickness, 
                       min_y=-width * 2, 
                       max_y=width * 2, 
                       field_function=None)

#     phi_w = get_weighting_potential(x, y, D=thickness, S=width, is_planar=True)
# 
#     # Plot weighting potential with colour and contour lines
#     levels = np.arange(0, 1, 5)
#     plt.imshow(phi_w, extent=[x.min(), x.max(), y.max(), y.min()], origin='upper', cmap=cm.get_cmap('Blues'))
#     plt.contour(x, y, phi_w, 15, colors='black')
# 
#     # Plot electrode
#     plt.plot([-width / 2., width / 2.], [0, 0], linewidth=5, color='red')
# 
#     # Plot weighting field directions
#     y, x = np.mgrid[0:thickness:20j, -width * 2:width * 2:20j]
#     E_x, E_y = get_weighting_field(x, y, D=thickness, S=width, is_planar=True)
# 
#     plt.quiver(x, y, E_x / np.sqrt(E_x ** 2 + E_y ** 2), E_y / np.sqrt(E_x ** 2 + E_y ** 2), pivot='mid', color='gray', scale=30.)
# 
#     plt.title('Weighting potential and weighting field direction (planar sensor)')
#     plt.xlabel('Position [um]')
#     plt.ylabel('Depth [um]')
#     plt.gca().set_aspect(1. / plt.gca().get_data_ratio() / 1.618)
#     plt.savefig('WeightingField_planar.pdf', layout='tight')
#     plt.show()