import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.patches import Rectangle
from matplotlib import cm
import matplotlib.animation as animation
import logging
_LOGGER = logging.getLogger(__name__)

from scarce import geometry


def get_mesh_plot(fig, mesh, values=None, invert_y_axis=True):
    ax = fig.add_subplot(111)

    vertexIDs = mesh._orderedCellVertexIDs
    vertexCoords = mesh.vertexCoords
    xCoords = np.take(vertexCoords[0], vertexIDs)
    yCoords = np.take(vertexCoords[1], vertexIDs)

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
                       depletion_function=None,
                       mesh=None,
                       title=None):
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
                           depletion_function=depletion_function,
                           mesh=mesh,
                           title=title)
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
                           depletion_function=None,
                           mesh=None,
                           title=None):
    _LOGGER.info('Create planar sensor plot')
    ax = fig.add_subplot(111)

    min_x, max_x = -width * float(n_pixel) / 2., width * float(n_pixel) / 2.
    min_y, max_y = 0., thickness

    # Define plot space
    x = np.linspace(min_x, max_x, 1000)
    y = np.linspace(min_y, max_y, 1000)
#     y = np.square(np.linspace(np.sqrt(min_y), np.sqrt(max_y), 1000)) TODO: make this work with streamplot

    # Create x,y plot grid
    xx, yy = np.meshgrid(x, y, sparse=True)

    def get_depletion_mask(depletion_function, x, y):
        ''' Returns true for all points outside of the depletion zone.
            If x, y in a sparse array representation are handled correctly.
        '''

        if x.shape != y.shape:  # Shapes do not fit
            if x.shape != y.T.shape:   # Sparse meshgrid assumption
                raise RuntimeError('The point representation in x,y in neither a grid nor a sparse grid.')
            x_dense, y_dense = np.meshgrid(x[0, :], y[:, 0], sparse=False)
        else:
            x_dense, y_dense = x, y

        mask = np.zeros_like(x_dense, dtype=np.bool)
        mask[y_dense > depletion_function(x_dense)] = True

        return mask

    if depletion_function:
        ax.plot(x, depletion_function(x), '-', color='black', linewidth=4)
        ax.plot(x, depletion_function(x), '--', color='blue', linewidth=4, label='Depletion')

    if potential_function:
        # Plot Potential
        phi = potential_function(xx, yy)

        if depletion_function:  # Mask potential in not depleted region, otherwise contour plot goes crazy
            phi = np.ma.masked_array(phi, mask=get_depletion_mask(depletion_function, xx, yy))

        # BUG in matplotlib: aspect to be set to equal, otherwise contour plot wrong aspect ratio
        # http://stackoverflow.com/questions/28857673/wrong-aspect-ratio-for-contour-plot-with-python-matplotlib
        ax.set_aspect('equal')
        ax.contour(x, y, phi, 10, colors='black', vmin=V_backplane, vmax=V_readout)
        cmesh = ax.pcolormesh(x - np.diff(x)[0] / 2., y - np.diff(y)[0] / 2., phi, cmap=cm.get_cmap('Blues'), vmin=V_backplane, vmax=V_readout)
        cax = fig.add_axes([ax.get_position().xmax, 0.1, 0.05, ax.get_position().ymax - ax.get_position().ymin])
        fig.colorbar(cmesh, cax=cax, orientation='vertical')

    # Plot E-Field
    if field_function:
        E_x, E_y = field_function(xx, yy)
        if depletion_function:  # Do not plot field in not depleted region, since it is 0
            E_x = np.ma.masked_array(E_x, mask=get_depletion_mask(depletion_function, xx, yy))
            E_y = np.ma.masked_array(E_y, mask=get_depletion_mask(depletion_function, xx, yy))

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
        ax.plot([pixel_position - width / 2, pixel_position - width / 2], [min_y, max_y], '--', color='black', linewidth=4)
    ax.plot([pixel_position + width / 2, pixel_position + width / 2], [min_y, max_y], '--', color='black', linewidth=4)  # Last pixel end

    ax.set_ylim((- 0.02 * (ax.get_ylim()[1] - ax.get_ylim()[0]), 1.02 * max_y))
    ax.set_xlabel('Position x/y [um]', fontsize=22)
    ax.set_ylabel('Position z [um]', fontsize=22)
    if title:
        ax.set_title(title, fontsize=22)
    ax.invert_yaxis()


def plot_3D_sensor(width_x, width_y,
                   radius, nD,
                   n_pixel_x, n_pixel_y,
                   V_bias, V_readout,
                   potential_function=None,
                   field_function=None,
                   mesh=None,
                   title=None):
    fig = plt.figure()
    get_3D_sensor_plot(fig,
                       width_x, width_y,
                       radius, nD,
                       n_pixel_x, n_pixel_y,
                       V_bias, V_readout,
                       potential_function=potential_function,
                       field_function=field_function,
                       mesh=mesh,
                       title=title)
    plt.show()


def get_3D_sensor_plot(fig,
                       width_x, width_y,
                       radius, nD,
                       n_pixel_x, n_pixel_y,
                       V_bias, V_readout,
                       potential_function=None,
                       field_function=None,
                       mesh=None,
                       title=None):

    ax = fig.add_subplot(111)

    desc = geometry.SensorDescription3D(width_x, width_y, n_pixel_x, n_pixel_y, radius, nD)

    min_x, max_x, min_y, max_y = desc.get_array_corners()

    # Define plot space
    x = np.linspace(min_x, max_x, max(width_x * n_pixel_x, width_y * n_pixel_y))  # >= 0.5 um resolution
    y = np.linspace(min_y, max_y, max(width_x * n_pixel_x, width_y * n_pixel_y))  # >= 0.5 um resolution

    # Create x,y plot grid
    xx, yy = np.meshgrid(x, y, sparse=True)

    # BUG in matplotlib: aspect to be set to equal, otherwise contour plot wrong aspect ratio
    # http://stackoverflow.com/questions/28857673/wrong-aspect-ratio-for-contour-plot-with-python-matplotlib
    ax.set_aspect('equal')

    def get_column_mask(x, y):
        ''' Returns true for points inside a column.
            X, y as a sparse array representation are handled correctly.
        '''

        if x.shape != y.shape:  # Shapes do not fit
            if x.shape != y.T.shape:   # Sparse meshgrid assumption
                raise RuntimeError('The point representation in x,y in neither a grid nor a sparse grid.')
            x_dense, y_dense = np.meshgrid(x[0, :], y[:, 0], sparse=False)
        else:
            x_dense, y_dense = x, y
        # Reduce radius to prevent aliasing at round column edges
        return desc.position_in_column(x_dense, y_dense)

    if potential_function:
        # Plot Potential
        phi = potential_function(xx, yy)

        # Mask potential in columns, otherwise contour plot goes crazy
        phi_masked = np.ma.masked_array(phi, mask=get_column_mask(xx, yy))

        ax.contour(x, y, phi_masked, 10, colors='black')
        cmesh = ax.pcolormesh(x - np.diff(x)[0] / 2., y - np.diff(y)[0] / 2., phi, cmap=cm.get_cmap('Blues'), vmin=V_bias, vmax=V_readout)
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

    # Plot pixel bondaries in x
    for pos_x in desc.get_pixel_x_offsets():
        ax.plot([pos_x - width_x / 2., pos_x - width_x / 2.], [min_y, max_y], '--', color='black', linewidth=4)
    ax.plot([pos_x + width_x / 2., pos_x + width_x / 2.], [min_y, max_y], '--', color='black', linewidth=4)  # Last pixel end

    # Plot pixel bondaries in y
    for pos_y in desc.get_pixel_y_offsets():
        ax.plot([min_x, max_x], [pos_y - width_y / 2., pos_y - width_y / 2.], '--', color='black', linewidth=4)
    ax.plot([min_x, max_x], [pos_y + width_y / 2., pos_y + width_y / 2.], '--', color='black', linewidth=4)  # Last pixel end

    # Plot readout pillars
    for pos_x, pos_y in desc.get_ro_col_offsets():
        ax.add_patch(plt.Circle((pos_x, pos_y), radius, color="darkred", linewidth=0))

    # Plot full bias pillars
    for pos_x, pos_y in desc.get_center_bias_col_offsets():
        ax.add_patch(plt.Circle((pos_x, pos_y), radius, color="darkblue", linewidth=0))

    # Plot side bias pillars
    for pos_x, pos_y in desc.get_side_bias_col_offsets():
        ax.add_patch(plt.Circle((pos_x, pos_y), radius, color="darkblue", linewidth=0))

    # Plot edge bias pillars
    for pos_x, pos_y in desc.get_edge_bias_col_offsets():
        ax.add_patch(plt.Circle((pos_x, pos_y), radius, color="darkblue", linewidth=0))

    ax.set_xlim((1.05 * min_x, 1.05 * max_x))
    ax.set_ylim((1.05 * min_y, 1.05 * max_y))
    ax.set_xlabel('Position x [um]', fontsize=22)
    ax.set_ylabel('Position y [um]', fontsize=22)
    if title:
        ax.set_title(title, fontsize=22)


def animate_drift_diffusion(fig, pe, ph, dt):
    if not fig.get_axes():
        raise RuntimeError('Function has to be called with sensor plot figure')
    else:
        ax = fig.get_axes()[0]

    electrons, = ax.plot([], [], '.', label='Electrons')
    holes, = ax.plot([], [], '.', label='Holes')
    time_template = '%.1f ns'
    time_text = ax.text(0.01, 0.97, '', transform=ax.transAxes)

    def init():
        electrons.set_data([], [])
        holes.set_data([], [])
        time_text.set_text('')
        return electrons, holes, time_text

    def animate(i):
        electrons.set_data(pe[i, 0, :], pe[i, 1, :])
        holes.set_data(ph[i, 0, :], ph[i, 1, :])
        time_text.set_text(time_template % (i * dt))
        return electrons, holes, time_text

    return init, animate


if __name__ == '__main__':
    import meshio as mio
    import fipy
    width_x, width_y = 250., 50.
    n_pixel_x, n_pixel_y = 3, 3
    radius, nD = 6., 2
    resolution = 10

    mesh = geometry.mesh_3D_sensor(width_x=width_x,
                                   width_y=width_y,
                                   n_pixel_x=n_pixel_x,
                                   n_pixel_y=n_pixel_y,
                                   radius=radius,
                                   nD=nD,
                                   resolution=resolution)

    plot_3D_sensor(width_x=width_x,
                   width_y=width_y,
                   radius=radius,
                   nD=nD,
                   n_pixel_x=n_pixel_x,
                   n_pixel_y=n_pixel_y,
                   V_bias=0,
                   V_readout=1,
                   potential_function=None,
                   field_function=None,
                   mesh=mesh,
                   title=None)


#     from scarce import fields
#
#     thickness = 200
#     width = 50
#     pitch = 50
#     n_pixel = 5
#     thickness = 200
#     resolution = 50
#     V_backplane, V_readout = -1, 0
#
#     def potential_function(x, y):
#         return fields.get_weighting_potential_analytic(x, y, D=thickness, S=width, is_planar=True)
#
#     def field_function(x, y):
#         return fields.get_weighting_field_analytic(x, y, D=thickness, S=width, is_planar=True)
#
#     plot_planar_sensor(width=width,
#                        pitch=width,
#                        thickness=thickness,
#                        n_pixel=n_pixel,
#                        V_backplane=0,
#                        V_readout=1,
#                        potential_function=potential_function,
#                        field_function=field_function,
#                        mesh=None)
