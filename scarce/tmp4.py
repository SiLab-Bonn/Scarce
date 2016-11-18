import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np

from scarce import geometry, fields

from fipy import GmshImporter2D


# Influences how correct the field for the center pixel(s) is
# due to more far away infinite boundary condition
n_pixel = 11

width = 50.
thickness = 200.

# Analytical solution only existing for pixel width = readout pitch (100% fill factor)
pitch = width

print 'width, thickness', width, thickness
geometry.mesh_planar_sensor(
    n_pixel=n_pixel,
    width=width,
    thickness=thickness,
    resolution=1000. * np.sqrt(width / 50.) * np.sqrt(50. / thickness),
    filename='planar_mesh_tmp_2.msh')
mesh = GmshImporter2D('planar_mesh_tmp_2.msh')

potential = fields.calculate_planar_sensor_w_potential(mesh=mesh,
                                                       width=width,
                                                       pitch=pitch,
                                                       n_pixel=n_pixel,
                                                       thickness=thickness)

potential_function = geometry.interpolate_potential(potential, smoothing=0.1)


def potential_analytic(x, y):
    return fields.get_weighting_potential_analytic(x, y, D=thickness, S=width, is_planar=True)

min_x, max_x = -width * float(n_pixel), width * float(n_pixel)
min_y, max_y = 0., thickness

nx, ny = 100 * n_pixel, 500
x = np.linspace(min_x, max_x, nx)
y = np.linspace(min_y, max_y, ny)

# Create x,y plot grid
xx, yy = np.meshgrid(x, y, sparse=True)

# Evaluate potential on a grid
pot_analytic = potential_analytic

E_w_y, E_w_x = np.gradient(-pot_analytic(xx, yy), np.diff(y)[0], np.diff(x)[0])
E_x_i, E_y_i = geometry.calculate_field(x, y, potential_function)

description = fields.Description(potential, 
                                 min_x=-width * float(n_pixel), 
                                 max_x=width * float(n_pixel), 
                                 min_y=0, 
                                 max_y=thickness,
                                 nx=200 * n_pixel,
                                 ny=200,
                                 smoothing=0.1
                                 )


for i in range(0, 100, 30):
#     plt.plot(y, pot_analytic(x[nx / 2 + i], y), label='Pot analytic')
#     plt.plot(y, potential_function(x[nx / 2 + i], y)[0, :],'.-', label='Pot numeric')
#     plt.legend(loc=1)

#     plt.twinx(plt.gca())
    E_x_test, E_y_test = description.get_field(x[nx / 2 + i], y)
    
    print E_y_test.shape

    plt.plot(y, E_w_y.T[nx / 2 + i, :], label='E analytic')
    plt.plot(y, E_y_i(x[nx / 2 + i], y)[0,:], label='E numeric')
    
    plt.plot(y, E_y_test[0,:], label='E numeric_2')

plt.legend(loc=0)
plt.show()
