import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Rectangle
from matplotlib import cm
import logging

import matplotlib
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'


from scarce import geometry, tools

_LOGGER = logging.getLogger(__name__)


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
    collection.set_edgecolor((0.2, 0.2, 0.2, 0.5))
    if invert_y_axis:
        ax.invert_yaxis()
    ax.add_collection(collection)

    if values:
        pass
        #     rgba = cm.cmap(plt.colors.norm(Z))
#         rgba = cm.get_cmap('coolwarm')(colors.norm(Z))
#         collection.set_facecolors(rgba)
#         collection.set_edgecolors(rgba)

#    ax.plot()


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
                       pot_func=None,
                       field_func=None,
                       depl_func=None,
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
                           pot_func=pot_func,
                           field_func=field_func,
                           depl_func=depl_func,
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
                           pot_func=None,
                           field_func=None,
                           depl_func=None,
                           mesh=None,
                           title=None):
    _LOGGER.info('Create planar sensor plot')
    ax = fig.add_subplot(111)

    min_x, max_x = -width * float(n_pixel) / 2., width * float(n_pixel) / 2.
    min_y, max_y = 0., thickness

    # Define plot space
    x = np.linspace(min_x, max_x, 1000)
    y = np.linspace(min_y, max_y, 1000)

#     TODO: make this work with streamplot
#     y = np.square(np.linspace(np.sqrt(min_y), np.sqrt(max_y), 1000))

    # Create x,y plot grid
    xx, yy = np.meshgrid(x, y, sparse=True)

    def get_depl_mask(depl_func, x, y):
        ''' Returns true for all points outside of the depletion zone.
            If x, y in a sparse array representation are handled correctly.
        '''

        if x.shape != y.shape:  # Shapes do not fit
            if x.shape != y.T.shape:   # Sparse meshgrid assumption
                raise RuntimeError('The point representation in x,y in neither'
                                   ' a grid nor a sparse grid.')
            x_dense, y_dense = np.meshgrid(x[0, :], y[:, 0], sparse=False)
        else:
            x_dense, y_dense = x, y

        mask = np.zeros_like(x_dense, dtype=np.bool)
        mask[y_dense > depl_func(x_dense)] = True

        return mask

    if depl_func:
        ax.plot(x, depl_func(x), '-', color='black', linewidth=4)
        ax.plot(x, depl_func(x), '--', color='blue', linewidth=4,
                label='Depletion')

    if pot_func:
        # Plot Potential
        phi = pot_func(xx, yy)

        # Mask pot in not depleted region, otherwise contour plot
        # goes crazy
        if depl_func:
            phi = np.ma.masked_array(
                phi, mask=get_depl_mask(depl_func, xx, yy))

        # BUG in matplotlib: aspect to be set to equal, otherwise contour plot
        # wrong aspect ratio
        # http://stackoverflow.com/questions/28857673/wrong-aspect-ratio-for-contour-plot-with-python-matplotlib
        ax.set_aspect('equal')
        ax.contour(x, y, phi, 10, colors='black', vmin=V_backplane,
                   vmax=V_readout)
        cmesh = ax.pcolormesh(x - np.diff(x)[0] / 2., y - np.diff(y)[0] / 2.,
                              phi, cmap=cm.get_cmap('Blues'), vmin=V_backplane,
                              vmax=V_readout, rasterized=True)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(cmesh, cax=cax, orientation='vertical')

    # Plot E-Field
    if field_func:
        E_x, E_y = field_func(xx, yy)
        # Do not plot field in not depleted region, since it is 0
        if depl_func:
            E_x = np.ma.masked_array(E_x, mask=get_depl_mask(depl_func,
                                                             xx, yy))
            E_y = np.ma.masked_array(E_y, mask=get_depl_mask(depl_func,
                                                             xx, yy))

        ax.streamplot(x, y, E_x, E_y, density=1.0, color='gray',
                      arrowstyle='-')
    elif pot_func:
        # Check for constant gradient
        assert np.allclose(np.gradient(x), np.gradient(x)[0])
        assert np.allclose(np.gradient(y), np.gradient(y)[0])
        E_y, E_x = np.gradient(-phi, np.gradient(y)[0], np.gradient(x)[0])
        ax.streamplot(x, y, E_x, E_y, density=1.0, color='gray',
                      arrowstyle='-')

    if mesh:
        get_mesh_plot(fig, mesh, invert_y_axis=False)

    # Plot backside
    ax.add_patch(Rectangle((min_x, max_y), (max_x - min_x), max_y,
                           color="darkblue", linewidth=2))

    # Plot pixel(s)
    ylim = ax.get_ylim()[0]
    for pixel in range(n_pixel):
        pixel_position = width * (pixel + 1. / 2.) - width * n_pixel / 2.
        ax.add_patch(Rectangle((pixel_position - pitch / 2, ylim),
                               pitch, 0, color="darkred", linewidth=5))
        ax.plot([pixel_position - width / 2, pixel_position - width / 2],
                [min_y, max_y], '--', color='black', linewidth=4)
    ax.plot([pixel_position + width / 2, pixel_position + width / 2],
            [min_y, max_y], '--', color='black', linewidth=4)  # Last pixel end

    ax.set_ylim((- 0.02 * (ax.get_ylim()[1] - ax.get_ylim()[0]), 1.02 * max_y))
    ax.set_xlabel('Position x/y [um]', fontsize=18)
    ax.set_ylabel('Position z [um]', fontsize=18)
    if title:
        ax.set_title(title, fontsize=18)
    ax.invert_yaxis()


def plot_3D_sensor(width_x, width_y,
                   radius, nD,
                   n_pixel_x, n_pixel_y,
                   V_bias, V_readout,
                   pot_func=None,
                   field_func=None,
                   mesh=None,
                   title=None):
    fig = plt.figure()
    get_3D_sensor_plot(fig,
                       width_x, width_y,
                       radius, nD,
                       n_pixel_x, n_pixel_y,
                       V_bias, V_readout,
                       pot_func=pot_func,
                       field_func=field_func,
                       mesh=mesh,
                       title=title)
    plt.show()


def get_3D_sensor_plot(fig,
                       width_x, width_y,
                       radius, nD,
                       n_pixel_x, n_pixel_y,
                       V_bias=None, V_readout=None,
                       pot_func=None,
                       field_func=None,
                       mesh=None,
                       title=None):

    ax = fig.add_subplot(111)

    desc = geometry.SensorDescription3D(width_x, width_y, n_pixel_x, n_pixel_y,
                                        radius, nD)

    min_x, max_x, min_y, max_y = desc.get_array_corners()

    # Define plot space with # >= 0.5 um resolution
    x = np.linspace(min_x, max_x, max(width_x * n_pixel_x,
                                      width_y * n_pixel_y))
    y = np.linspace(min_y, max_y, max(width_x * n_pixel_x,
                                      width_y * n_pixel_y))

    # Create x,y plot grid
    xx, yy = np.meshgrid(x, y, sparse=True)

    # BUG in matplotlib: aspect to be set to equal, otherwise contour plot
    # has wrong aspect ratio
    # http://stackoverflow.com/questions/28857673/wrong-aspect-ratio-for-contour-plot-with-python-matplotlib
    ax.set_aspect('equal')

    def get_column_mask(x, y):
        ''' Returns true for points inside a column.
            X, y as a sparse array representation are handled correctly.
        '''

        if x.shape != y.shape:  # Shapes do not fit
            if x.shape != y.T.shape:   # Sparse meshgrid assumption
                raise RuntimeError('The point representation in x,y in neither\
                 a grid nor a sparse grid.')
            x_dense, y_dense = np.meshgrid(x[0, :], y[:, 0], sparse=False)
        else:
            x_dense, y_dense = x, y
        # Reduce radius to prevent aliasing at round column edges
        return desc.position_in_column(x_dense, y_dense)

    if pot_func:
        # Plot Potential
        phi = pot_func(xx, yy)

        # Mask pot in columns, otherwise contour plot goes crazy
        phi_masked = np.ma.masked_array(phi, mask=get_column_mask(xx, yy))

        if V_bias is None:
            V_bias = phi_masked.min()

        if V_readout is None:
            V_readout = phi_masked.max()

        ax.contour(x, y, phi_masked, 10, colors='black')
        cmesh = ax.pcolormesh(x - np.diff(x)[0] / 2., y - np.diff(y)[0] / 2.,
                              phi, cmap=cm.get_cmap('Blues'), vmin=V_bias,
                              vmax=V_readout, rasterized=True)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(cmesh, cax=cax, orientation='vertical')

    # Plot E-Field
    if field_func:
        E_x, E_y = field_func(xx, yy)
        ax.streamplot(x, y, E_x, E_y, density=1.0, color='gray',
                      arrowstyle='-')
    elif pot_func:  # Get field from pot differentiation
        assert np.allclose(np.gradient(x), np.gradient(x)[0])
        assert np.allclose(np.gradient(y), np.gradient(y)[0])
        E_y, E_x = np.gradient(-phi, np.gradient(y)[0], np.gradient(x)[0])
        ax.streamplot(x, y, E_x, E_y, density=1.0, color='gray',
                      arrowstyle='-')

    if mesh:
        get_mesh_plot(fig, mesh, invert_y_axis=False)

    # Plot pixel bondaries in x
    for pos_x in desc.get_pixel_x_offsets():
        ax.plot([pos_x - width_x / 2., pos_x - width_x / 2.], [min_y, max_y],
                '--', color='black', linewidth=4)
    # Last pixel end
    ax.plot([pos_x + width_x / 2., pos_x + width_x / 2.], [min_y, max_y],
            '--', color='black', linewidth=4)

    # Plot pixel bondaries in y
    for pos_y in desc.get_pixel_y_offsets():
        ax.plot([min_x, max_x], [pos_y - width_y / 2., pos_y - width_y / 2.],
                '--', color='black', linewidth=4)
    ax.plot([min_x, max_x], [pos_y + width_y / 2., pos_y + width_y / 2.],
            '--', color='black', linewidth=4)  # Last pixel end

    # Plot readout pillars
    for pos_x, pos_y in desc.get_ro_col_offsets():
        ax.add_patch(plt.Circle((pos_x, pos_y), radius, color="darkred",
                                linewidth=0, zorder=5))

    # Plot full bias pillars
    for pos_x, pos_y in desc.get_center_bias_col_offsets():
        ax.add_patch(plt.Circle((pos_x, pos_y), radius, color="darkblue",
                                linewidth=0, zorder=5))

    # Plot side bias pillars
    for pos_x, pos_y in desc.get_side_bias_col_offsets():
        ax.add_patch(plt.Circle((pos_x, pos_y), radius, color="darkblue",
                                linewidth=0, zorder=5))

    # Plot edge bias pillars
    for pos_x, pos_y in desc.get_edge_bias_col_offsets():
        ax.add_patch(plt.Circle((pos_x, pos_y), radius, color="darkblue",
                                linewidth=0, zorder=5))

    ax.set_xlim((1.05 * min_x, 1.05 * max_x))
    ax.set_ylim((1.05 * min_y, 1.05 * max_y))
    ax.set_xlabel('Position x [um]', fontsize=18)
    ax.set_ylabel('Position y [um]', fontsize=18)
    if title:
        ax.set_title(title, fontsize=18)


def animate_drift_diffusion(fig, T, pe, ph, dt, n_steps):
    if not fig.get_axes():
        raise RuntimeError('Function has to be called with sensor plot figure')
    else:
        ax = fig.get_axes()[0]

    # Interpolate to common time base
    pe_x = tools.time_data_interpolate(T, pe[:, 0],
                                       t=np.linspace(
        0, dt * n_steps, n_steps),
        fill_value=np.nan)
    pe_y = tools.time_data_interpolate(T, pe[:, 1],
                                       t=np.linspace(
        0, dt * n_steps, n_steps),
        fill_value=np.nan)
    ph_x = tools.time_data_interpolate(T, ph[:, 0],
                                       t=np.linspace(
        0, dt * n_steps, n_steps),
        fill_value=np.nan)
    ph_y = tools.time_data_interpolate(T, ph[:, 1],
                                       t=np.linspace(
        0, dt * n_steps, n_steps),
        fill_value=np.nan)

    electrons, = ax.plot([], [], 'o', label='Electrons')
    holes, = ax.plot([], [], '.', label='Holes')
    time_template = '%.1f ns'
    time_text = ax.text(0.01, 0.97, '', transform=ax.transAxes)

    def init():
        electrons.set_data([], [])
        holes.set_data([], [])
        time_text.set_text('')
        return electrons, holes, time_text

    def animate(i):
        electrons.set_data(pe_x[i], pe_y[i])
        holes.set_data(ph_x[i], ph_y[i])
        time_text.set_text(time_template % (i * dt))
        return electrons, holes, time_text

    return init, animate
