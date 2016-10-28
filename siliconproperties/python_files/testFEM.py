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

from siliconproperties.python_files.plot import plot_mesh
from matplotlib import colors, cm


def generate_3D_pixel(geom, x=250., y=50., radius=6., nD=2, resolution=10.):
    pillars = []

    # Create readout pillars
    for pillar in range(nD):
        position = x / nD * (pillar + 1. / 2.) - x / 2.
        circle = geom.add_circle(x0=[position, 0.0, 0.0],
                                 radius=radius,
                                 lcar=resolution / 2.,
                                 num_sections=4,
                                 # If compound==False, the section borders have to be points of the
                                 # discretization. If using a compound circle, they don't; gmsh can
                                 # choose by itself where to point the circle points.
                                 compound=False
                                 )
        pillars.append(geom.add_line_loop(circle))


def mesh_3D_sensor(x=250., y=50., radius=6., nD=2, resolution=10.):
    geom = pg.Geometry()

    pillars = []

    generate_3D_pixel(geom, x, y, radius, nD, resolution)

    # Create readout
    print -y / 2
    circle_points = []
    circle_points.append(geom.add_point([-x / 4, -y / 2, 0], lcar=resolution))
    circle_points.append(geom.add_point([0, y, 0], lcar=resolution))
    circle_points.append(geom.add_point([x / 4, -y / 2, 0], lcar=resolution))
#     circle_half = geom.add_circle_sector(circle_points)

    t = [circle_points[0] + circle_points[2]]
    #circle_half_ll = geom.add_compound_line(circle_points)

    points = []

    points_xyz = [
        #[-x/4, -y/4, 0],
        #[x/4, -y/4, 0],
        #[x/4, -y/2, 0],
        [x / 2, -y / 2, 0],
        [x / 2, y / 2, 0],
        [-x / 2, y / 2, 0],
        [-x / 2, -y / 2, 0],
        #[-x/4, -y/2, 0]
    ]
    for point in points_xyz:
        points.append(geom.add_point(point, lcar=resolution))

#     points = t + points

    # Create lines
    lines = [geom.add_line(points[i], points[i + 1])
             for i in range(len(points) - 1)]
    lines.append(geom.add_line(points[-1], points[0]))

    line_loop = geom.add_line_loop(lines)
    geom.add_plane_surface([line_loop])

    return geom


def mesh_planar_sensor(x, thickness, resolution=1.):
    geom = pg.Geometry()
    resolution_x = x / resolution

    points_xyz = [
        [x / 2, -thickness / 2, 0],
        [x / 2, thickness / 2, 0],
        [-x / 2, thickness / 2, 0],
        [-x / 2, -thickness / 2, 0],
    ]

    points = []
    points.append(geom.add_point(points_xyz[0], lcar=resolution_x))
    points.append(geom.add_point(points_xyz[1], lcar=resolution_x))
    points.append(geom.add_point(points_xyz[2], lcar=resolution_x))
    points.append(geom.add_point(points_xyz[3], lcar=resolution_x))

    # Create lines
    lines = [geom.add_line(points[i], points[i + 1])
             for i in range(len(points) - 1)]
    lines.append(geom.add_line(points[-1], points[0]))

    line_loop = geom.add_line_loop(lines)
    geom.add_plane_surface([line_loop])

    # Add 1/x1.5 law for the mesh size
    raw_codes = ['lc = %f;' % (resolution_x / 4.),
                 'Field[1] = Attractor;',
                 'Field[1].EdgesList = {l2};'
                 'Field[1].NNodesByEdge = %d;' % resolution,
                 'Field[2] = MathEval;',
                 'Field[2].F = Sprintf(\"F1^3 + %g\", lc);',
                 'Background Field = 2;\n']

    geom.add_raw_code(raw_codes)
    return geom


def calculate_planar_sensor_potential(width, pitch, n_pixel, thickness,
                                      resolution, V_backplane, V_readout=0):
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

    electrodes = readout_plane
    bcs = [fipy.FixedValue(value=V_backplane, faces=backplane)]
    X, _ = mesh.getFaceCenters()
    for pixel in range(n_pixel):
        pixel_position = width * (pixel + 1. / 2.) - width * n_pixel / 2.
        bcs.append(fipy.FixedValue(value=V_readout,
                                   faces=electrodes &
                                   (X > pixel_position - pitch / 2.) &
                                   (X < pixel_position + pitch / 2.)))

    potential.equation.solve(var=potential, boundaryConditions=bcs)
    return potential


def interpolate_potential(potential, smoothing):
    x = np.array(potential.mesh.getFaceCenters()[0])
    y = np.array(potential.mesh.getFaceCenters()[1])
    z = np.array(potential.arithmeticFaceValue)
    return SmoothBivariateSpline(x, y, z, s=smoothing, kx=3, ky=3)

if __name__ == '__main__':
    width = 50
    pitch = 45
    n_pixel = 3
    thickness = 50
    resolution = 100.
    V_backplane, V_readout = -1, 0

    potential = calculate_planar_sensor_potential(width, pitch, n_pixel, thickness, resolution, V_backplane, V_readout)

    plot_mesh(potential.mesh)
 
    min_x, max_x = np.min(np.array(potential.mesh.getFaceCenters()[0])), np.max(np.array(potential.mesh.getFaceCenters()[0]))
    min_y, max_y = np.min(np.array(potential.mesh.getFaceCenters()[1])), np.max(np.array(potential.mesh.getFaceCenters()[1]))
 
    print 'Interpolate'
    fit = interpolate_potential(potential, smoothing=1)
 
    # Plot potential
    xnew = np.linspace(min_x, max_x, 1000)
    ynew = np.linspace(min_y, max_y, 1000)
    phi = fit(xnew, ynew).T
    plt.contour(xnew, ynew, phi, 15, colors='black')
    plt.pcolormesh(xnew, ynew, phi, cmap=cm.get_cmap('Blues'), vmin=V_backplane, vmax=V_readout)
    plt.colorbar()
    
#     print phi.shape
#     plt.clf()
#     plt.plot(ynew, phi.T[600])
#     plt.plot(ynew, phi.T[500])
#     plt.plot(ynew, phi.T[400])
#     plt.show()
#     raise
 
    # Plot E-Field
    xnew = np.linspace(min_x, max_x, 10)
    ynew = np.linspace(min_y, max_y, 10)
    z_new = fit(xnew, ynew).T
    xnew, ynew = np.meshgrid(xnew, ynew)
    E_y, E_x = np.gradient(z_new)
    E_x, E_y = -E_x, -E_y
    plt.streamplot(xnew, ynew, E_x, E_y, density=1.5, color='gray', arrowstyle='-')
    plt.savefig('FEM_solution.pdf')
    plt.show()
    

