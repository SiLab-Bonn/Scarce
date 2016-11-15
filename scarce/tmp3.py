''' gmsh has to be installed and put into an environment variable (e.g. PATH)
    in a way that the command gmsh from the terminal/console starts it.
'''

import pygmsh as pg
import numpy as np
import meshio as mio
import fipy
import matplotlib.pyplot as plt

from scipy import interpolate
from scipy.interpolate import SmoothBivariateSpline, LSQBivariateSpline, UnivariateSpline, RectBivariateSpline

from scarce.plot import plot_mesh
from matplotlib import colors, cm
from scarce import geometry


def mesh_planar_sensor(x, thickness, resolution=3.):
    geom = pg.Geometry()

    geom.add_rectangle(xmin=-x / 2.,
                       xmax=x / 2.,
                       ymin=-thickness / 2.,
                       ymax=thickness / 2.,
                       z=0,
                       lcar=resolution)

    return geom


def calculate_planar_sensor_potential(width, pitch, n_pixel, thickness, resolution, V_backplane, V_readout=0):
    points, cells = pg.generate_mesh(mesh_planar_sensor(x=width * n_pixel,
                                                        thickness=thickness,
                                                        resolution=resolution))

    mio.write('sensor.msh', points, cells)
    mesh = fipy.GmshImporter2D('sensor.msh')

    potential = fipy.CellVariable(mesh=mesh, name='potential', value=0.)
    permittivity = 1.
    potential.equation = (fipy.DiffusionTerm(coeff=permittivity) == 0.)

    # Calculate boundaries
    V_backplane = V_backplane
    backplane = mesh.getFacesBottom()

    V_readout = V_readout
    readout_plane = mesh.getFacesTop()

    #no_fill = width - pitch
    electrodes = readout_plane
    bcs = [fipy.FixedValue(value=V_backplane, faces=backplane)]
    X, _ = mesh.getFaceCenters()
    for pixel in range(n_pixel):
        pixel_position = width * (pixel + 1. / 2.) - width * n_pixel / 2.
        bcs.append(fipy.FixedValue(value=V_readout,
                                   faces=electrodes & (X > pixel_position - pitch / 2.) & (X < pixel_position + pitch / 2.)))

    potential.equation.solve(var=potential, boundaryConditions=bcs)
    return potential


def interpolate_potential(x, y, z, smoothing):
    return SmoothBivariateSpline(x, y, z, s=smoothing, kx=3, ky=3)

if __name__ == '__main__':
    width = 50
    pitch = 5
    n_pixel = 9
    thickness = 200
    resolution = 1.5
    V_backplane, V_readout = -1, 0

    potential = calculate_planar_sensor_potential(width, pitch, n_pixel, thickness, resolution, V_backplane, V_readout)

#     plot_mesh(potential.mesh)

    print 'Interpolate'

    min_x, max_x = np.min(np.array(potential.mesh.getFaceCenters()[0])), np.max(np.array(potential.mesh.getFaceCenters()[0]))
    min_y, max_y = np.min(np.array(potential.mesh.getFaceCenters()[1])), np.max(np.array(potential.mesh.getFaceCenters()[1]))

    print np.sum(np.square(potential.value))
#     raise
    x = np.array(potential.mesh.getFaceCenters()[0])
    y = np.array(potential.mesh.getFaceCenters()[1])
    z = np.array(potential.arithmeticFaceValue)
    
    potential_inter = geometry.interpolate_potential(potential)
    xtest = np.linspace(min_x, max_x, 1000)
    ytest = np.linspace(min_y, max_y, 1000)
    xxtest, yytest = np.meshgrid(xtest, ytest)
    zztest = potential_inter(xxtest, yytest)
    sel = np.isfinite(zztest.ravel())
    fit_test = interpolate_potential(xxtest.ravel()[sel], yytest.ravel()[sel], zztest.ravel()[sel], smoothing=1000000)
#     interpolate.RectBivariateSpline(xtest, ytest, zztest)
    print xtest, ytest, zztest
#     raise
    print '_____'
    print x, y, z
    raise
    fit = interpolate_potential(x, y, z, smoothing=100)

    # Plot potential
    xnew = np.linspace(min_x, max_x, 1000)
    ynew = np.linspace(min_y, max_y, 1000)
    phi = fit(xnew, ynew).T
    plt.contour(xnew, ynew, phi, 15, colors='black')
    plt.pcolormesh(xnew, ynew, phi, cmap=cm.get_cmap('Blues'))
    plt.colorbar()

    # Plot E-Field
    xnew = np.linspace(min_x, max_x, 10)
    ynew = np.linspace(min_y, max_y, 10)
    z_new = fit(xnew, ynew).T
    xnew, ynew = np.meshgrid(xnew, ynew)
    E_y, E_x = np.gradient(z_new)
    E_x, E_y = -E_x, -E_y
    plt.streamplot(xnew, ynew, E_x, E_y, density=1.5, color='gray', arrowstyle='-')
    plt.gca().set_aspect(1.)
    plt.show()
