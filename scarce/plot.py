import matplotlib.pyplot as plt
from fipy.tools import numerix
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib import colors, cm


def plot_mesh(mesh, values=None, invert_y_axis=True):

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
    if invert_y_axis:
        plt.gca().invert_yaxis()
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
                       V_backplane, 
                       V_readout, 
                       min_x, max_x, 
                       min_y, max_y,
                       e_field_function=None):

    # Define plot space
    x = np.linspace(min_x, max_x, 1000)
    y = np.linspace(min_y, max_y, 1000)

    # Interpolate potential to a x,y grid
    xx, yy = np.meshgrid(x, y, sparse=True)
    phi = potential_function(xx, yy)

    # Plot Potential
    
    # BUG in matplotlib: aspect to be set to equal, otherwise contour plot wrong aspect ratio
    # http://stackoverflow.com/questions/28857673/wrong-aspect-ratio-for-contour-plot-with-python-matplotlib
    plt.gca().set_aspect('equal')
    plt.contour(x, y, phi, 10, colors='black')
    plt.pcolormesh(x, y, phi, cmap=cm.get_cmap('Blues'), vmin=V_backplane, vmax=V_readout)
    plt.colorbar()
    
    # Plot E-Field
    if e_field_function:
        x = np.linspace(min_x, max_x, 20)
        y = np.linspace(min_y, max_y, 20)
        xx, yy = np.meshgrid(x, y, sparse=True)
        E_x, E_y = e_field_function(xx, yy)
    else:
        E_x, E_y  = np.gradient(-phi, np.diff(x)[0], np.diff(y)[0], axis=(1, 0))

    plt.streamplot(x, y, E_x, E_y, density=1.0, color='gray', arrowstyle='-')
   
    # Plot pixel background
    plt.gca().add_patch(plt.Rectangle((min_x, plt.ylim()[1]), (max_x - min_x), 0.05*(plt.ylim()[1] - plt.ylim()[0]), color="grey", linewidth=0))
    
    # Plot pixel(s)
    for pixel in range(n_pixel):
        pixel_position = width * (pixel + 1. / 2.) - width * n_pixel / 2.
        plt.gca().add_patch(plt.Rectangle((pixel_position - pitch / 2, plt.ylim()[0]), pitch, 0.05*(plt.ylim()[1] - plt.ylim()[0]), color="darkred", linewidth=0))
    
    plt.ylim((1.05 * min_y, 1.05 * max_y))
    plt.xlabel('Position x/y [um]', fontsize=22)
    plt.ylabel('Position z [um]', fontsize=22)
    plt.gca().invert_yaxis()
    plt.show()
    

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
    positions = ([- pitch_x/2., -pitch_y/2.],
                 [0, -pitch_y/2.],
                 [pitch_x/2., -pitch_y/2.],
                 [- pitch_x/2., pitch_y/2.],
                 [0, pitch_y/2.],
                 [pitch_x/2., pitch_y/2.])
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
    width = 40  # [um]

    def potential_function(x, y):
        return fields.get_weighting_potential(x, y, D=thickness, S=width, is_planar=True)
    
    def e_field_function(x, y):
        return fields.get_weighting_field(x, y, D=thickness, S=width, is_planar=True)
    
    e_field_function = None
    
    plot_planar_sensor(potential_function=potential_function, 
                       width=width, 
                       pitch=width, 
                       n_pixel=1, 
                       V_backplane=0, 
                       V_readout=1, 
                       min_x=-width * 2, 
                       max_x=width * 2, 
                       min_y=0, 
                       max_y=thickness,
                       e_field_function=e_field_function)
