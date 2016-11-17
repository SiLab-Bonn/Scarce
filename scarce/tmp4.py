import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np

from scarce import geometry, fields

from fipy import GmshImporter2D


def interpolate_nan(a):
    mask = np.isnan(a)
    a[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), a[~mask])
    return a

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
    resolution=600. * np.sqrt(width / 50.) * np.sqrt(50. / thickness),
    filename='planar_mesh_tmp_2.msh')
mesh = GmshImporter2D('planar_mesh_tmp_2.msh')

potential = fields.calculate_planar_sensor_w_potential(mesh=mesh,
                                                       width=width,
                                                       pitch=pitch,
                                                       n_pixel=n_pixel,
                                                       thickness=thickness)

potential_function = geometry.interpolate_potential(potential)


def potential_analytic(x, y):
    return fields.get_weighting_potential_analytic(x, y, D=thickness, S=width, is_planar=True)

min_x, max_x = -width * float(n_pixel), width * float(n_pixel)
min_y, max_y = 0., thickness

nx, ny = 200 * n_pixel, 500
x = np.linspace(min_x, max_x, nx)
y = np.linspace(min_y, max_y, ny)

# Create x,y plot grid
xx, yy = np.meshgrid(x, y, sparse=True)

# Evaluate potential on a grid
pot_analytic = potential_analytic
pot_numeric = potential_function(xx, yy)

pot_numeric = interpolate_nan(pot_numeric)

# Transforms from pot(y, x) to pot(x, y)
pot_numeric_smooth = interpolate.RectBivariateSpline(x, y, pot_numeric.T, s=0.1)

print np.any(~np.isfinite(pot_numeric_smooth(xx, yy)))
#     raise

E_w_y, E_w_x = np.gradient(-pot_analytic(xx, yy), np.diff(y)[0], np.diff(x)[0])
E_w_x_2, E_w_y_2 = np.gradient(-pot_numeric_smooth(xx, yy), np.diff(x)[0], np.diff(y)[0])
#
#     print E_w_y, E_w_x

i = 0
plt.plot(y, pot_analytic(xx, yy).T[nx / 2 + i, :], label='Pot analytic')
plt.plot(y, pot_numeric_smooth(xx, yy)[nx / 2 + i, :], label='Pot numeric smooth')
#     i = 40
#     plt.plot(y, pot_analytic(xx, yy).T[nx / 2 + i, :], label='Pot analytic')
#     plt.plot(y, pot_numeric_smooth(xx, yy)[nx / 2 + i, :], label='Pot numeric smooth')

plt.twinx(plt.gca())

plt.plot(y, E_w_y.T[nx / 2 + i, :], label='Analytic diff')
plt.plot(y, E_w_y_2[nx / 2 + i, :], label='Numeric diff')
#     plt.ylim((0, plt.ylim()[1] / 2.))
#     plt.plot(y, pot_numeric.T[nx / 2 + i, :], label='Numeric')
#     plt.plot(y, tt(np.ones_like(y) * x[nx / 2 + i], y[:])[0], label='Numeric_inter')
plt.legend(loc=0)
plt.show()

# for i in [-45, -30, -15, -10, 0, 10, 15, 30, 45]:  # Check only at center pixel, edge pixel are not interessting
#         sel = pot_analytic.T[nx / 2 + i, :] > 0.01
# Check with very tiny and tuned error allowance
#         self.assertTrue(np.allclose(pot_analytic.T[nx / 2 + i, sel], pot_numeric.T[nx / 2 + i, sel], rtol=0.01, atol=0.005))

raise
width = 50.
# Analytical solution only existing for pixel width = readout pitch (100 % fill factor)
pitch = width
thickness = 200.
n_pixel = 11

geometry.mesh_planar_sensor(
    n_pixel=n_pixel,
    width=width,
    thickness=thickness,
    resolution=200,
    filename='planar_mesh_tmp_2.msh')
mesh = GmshImporter2D('planar_mesh_tmp_2.msh')

potential = fields.calculate_planar_sensor_w_potential(mesh=mesh,
                                                       width=width,
                                                       pitch=pitch,
                                                       n_pixel=n_pixel,
                                                       thickness=thickness)

field_function = geometry.calculate_field(potential)


def field_analytic(x, y):
    return fields.get_weighting_field_analytic(x, y, D=thickness, S=width, is_planar=True)

min_x, max_x = -width * float(n_pixel), width * float(n_pixel)
min_y, max_y = 0., thickness

nx, ny = 100, 100
x = np.linspace(min_x, max_x, nx)
y = np.linspace(min_y, max_y, ny)

# Create x,y plot grid
xx, yy = np.meshgrid(x, y, sparse=False)

# Evaluate potential on a grid
f_analytic = field_analytic(xx, yy)
f_numeric = field_function(xx, yy)

print x.shape, y.shape, f_analytic[1].shape
print len(xx.ravel()), len(yy.ravel()), len(f_analytic[1].ravel())


f_y = f_numeric[1].T
f_y_a = f_analytic[1].T

tt = interpolate.RectBivariateSpline(x, y, f_y_a)

f_y = interpolate_nan(f_y)
print f_y
assert np.count_nonzero(~np.isfinite(f_y)) == 0
#     f_y[~np.isfinite(f_y)] = 0.

#     sel = np.where(np.isfinite(f_y))
#     print sel[0], np.any(np.diff(x[sel[0]] <0.))
#     print sel[1], np.any(np.diff(y[sel[1]] <0.))

ttt = interpolate.RectBivariateSpline(x[1:-1], y[1:-1], f_y[1:-1, 1:-1], s=0.01)

#     for i in range(1000):
#         ttt(x[nx / 2], y[nx / 2])[0]
#         print i

#     tttt = interpolate.SmoothBivariateSpline(xx.ravel(), yy.ravel(), f_y.ravel())
#     raise
#     print np.count_nonzero(~np.isfinite(f_analytic[1])), np.count_nonzero(~np.isfinite(f_numeric[1]))
#     ttt = interpolate.SmoothBivariateSpline(xx.ravel(), yy.ravel(), f_analytic[1].ravel())

#     ttt = interpolate.bisplrep(xx, yy, f_analytic[1])

i, j = 0, 50
print f_y_a[nx / 2, :]
print tt(np.ones_like(y) * x[nx / 2], y[:])[0]
print ttt(np.ones_like(y) * x[nx / 2], y[:])[0]

#     print np.allclose(f_analytic[1].T[nx / 2, :], tt(np.ones_like(y) * x[nx / 2], y[:])[nx / 2])
#     raise
#     print ttt
#     print ttt(x[i], y[j])

#     print tt(0, y[1:])[0]
#
i = 0
plt.plot(y[:], f_y_a[nx / 2 + i, :], label='Analytic')
plt.plot(y[:], f_y[nx / 2 + i, :], label='Numeric')
plt.plot(y[:], ttt(np.ones_like(y) * x[nx / 2 + i], y[:])[0], label='Numeric interpolate')
plt.plot(y[:], tt(np.ones_like(y) * x[nx / 2 + i], y[:])[0], label='Analytic interpolate')
plt.legend(loc=0)
plt.show()
raise
for i in [-45, -30, -15, -10, 0, 10, 15, 30, 45]:  # Check only at center pixel, edge pixel are not interessting
    sel = f_analytic[0].T[nx / 2 + i, :] > 0.001
    assert(np.allclose(f_analytic[0].T[nx / 2 + i, sel], f_numeric[0].T[nx / 2 + i, sel], rtol=0.01, atol=0.01))
    sel = f_analytic[1].T[nx / 2 + i, :] > 0.001
    assert(np.allclose(f_analytic[1].T[nx / 2 + i, sel], f_numeric[1].T[nx / 2 + i, sel], rtol=0.01, atol=0.03))
