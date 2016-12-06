''' Example that creates a planar silicon sensor with a given geometry (thickness, number of pixels, pitch, width).
    Calculates the electrical potential and fields. For comparison also the analytical result of a planar sensor with
    100% fill factor (width = pitch) is created.

    .. WARNING::
       The calculation of the depletion region is simplified. If the depletion is not at a contant y position
       in the sensor (e.g. for pixels with small fill factor) it deviates from the correct solution.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scarce import fields, plot, geometry, silicon, tools, solver


def induced_current():
    # Number of pixels influences how correct the field for the
    # center pixel(s) is due to more far away infinite boundary condition
    n_pixel = 9

    # Geometry of one pixel
    width = 50.
    thickness = 300.
    pitch = 45.

    n_eff = 5e12  # n_eff [cm^-3]
    temperature = 300

    # Potentials
    V_bias = -160.
    V_readout = 0.
    V_bi = -silicon.get_diffusion_potential(n_eff, temperature)

    # Create mesh of the sensor and stores the result
    # The created file can be viewed with any mesh viewer (e.g. gmsh)
    mesh = geometry.mesh_planar_sensor(
        n_pixel=n_pixel,
        width=width,
        thickness=thickness,
        resolution=200.,
        filename='planar_mesh_example.msh')

    # Numerically solve the laplace equation on the mesh
    potential = fields.calculate_planar_sensor_potential(mesh=mesh,
                                                         width=width,
                                                         pitch=pitch,
                                                         n_pixel=n_pixel,
                                                         thickness=thickness,
                                                         n_eff=n_eff,
                                                         V_bias=V_bias,
                                                         V_readout=V_readout,
                                                         V_bi=V_bi)

    w_potential = fields.calculate_planar_sensor_w_potential(mesh=mesh,
                                                             width=width,
                                                             pitch=pitch,
                                                             n_pixel=n_pixel,
                                                             thickness=thickness)

    min_x = float(mesh.getFaceCenters()[0, :].min())
    max_x = float(mesh.getFaceCenters()[0, :].max())
    pot_descr = fields.Description(potential,
                                   min_x=min_x,
                                   max_x=max_x,
                                   min_y=0,
                                   max_y=thickness)
    pot_w_descr = fields.Description(w_potential,
                                     min_x=min_x,
                                     max_x=max_x,
                                     min_y=0,
                                     max_y=thickness)

    pot_descr.get_field(0, 0)
    pot_w_descr.get_field(0, 0)

    tools.save(pot_descr, 'planar_sensor_pot.sc')
    tools.save(pot_w_descr, 'planar_sensor_pot_w.sc')
#     raise
    pot_descr = tools.load('planar_sensor_pot.sc')
    pot_w_descr = tools.load('planar_sensor_pot_w.sc')

    # Get analytical result
    def potential_analytic(x, y):
        return fields.get_weighting_potential_analytic(x, y, D=thickness, S=width, is_planar=True)

# Plot analytical / numerical result with depletion region in 1D
#     y = np.linspace(0, thickness, 100)
#     x = np.zeros_like(y)
#     plt.plot(y, pot_descr.get_potential(x, y), label='Potential, numerical', linewidth=2)
#     pot_masked = np.ma.masked_array(pot_descr.get_potential(x, y), mask=pot_descr.get_depletion_mask(x, y), linewidth=2)
#     plt.plot(y, pot_masked, label='Potential, numerical, depleted', linewidth=2)
#     plt.plot([pot_descr.get_depletion(x[50]), pot_descr.get_depletion(x[50])], plt.ylim(), label='Depletion, numerical ', linewidth=2)
#     plt.plot(y, fields.get_potential_planar_analytic_1D(y, V_bias=V_bias + V_bi, V_readout=V_readout, n_eff=n_eff, D=thickness), '--', label='Potential, analytical', linewidth=2)
#     plt.plot([silicon.get_depletion_depth(np.abs(V_bias), n_eff / 1e12, temperature), silicon.get_depletion_depth(np.abs(V_bias), n_eff / 1e12, temperature)], plt.ylim(), '--', label='Depletion, analytical', linewidth=2)
#     plt.legend(loc=0)
#     plt.title('Potential in a not fully depleted planar sensor')
#     plt.xlabel('Position [um]')
#     plt.xlabel('Potential [V]')
#     plt.grid()
#     plt.show()

# Plot numerical result in 2D
#     plot.plot_planar_sensor(width=width,
#                             pitch=pitch,
#                             thickness=thickness,
#                             n_pixel=n_pixel,
#                             V_backplane=0.,
#                             V_readout=1.,
#                             potential_function=pot_w_descr.get_potential,
#                             field_function=pot_w_descr.get_field,
#                             title='Planar sensor weighting potential')
#     plot.plot_planar_sensor(width=width,
#                             pitch=pitch,
#                             thickness=thickness,
#                             n_pixel=n_pixel,
# V_backplane=V_bias,  # Weighting field = 0 at backplane
# V_readout=V_readout,  # Weighting field = 1 at readout
#                             potential_function=pot_descr.get_potential,
#                             field_function=pot_descr.get_field,
#                             depletion_function=pot_descr.get_depletion,
#                             title='Planar sensor potential')

    xx, yy = np.meshgrid(np.linspace(0, width, 1.), np.linspace(thickness / 10., thickness, 1.), sparse=False)
    p0 = np.array([xx.ravel(), yy.ravel()])
    q0 = np.ones_like(p0)

    dd = solver.DriftDiffusionSolver(pot_descr, pot_w_descr, T=temperature)

    dt = 0.01
    n_steps = 800

    traj_e, traj_h, I_ind_e, I_ind_h = dd.solve(p0, q0, dt, n_steps)
#     traj_coarse = dd.solve(p0, dt * 10., n_steps / 10)

    plt.plot(np.arange(n_steps) * dt, I_ind_e[:, 0] / n_steps)
    plt.plot(np.arange(n_steps) * dt, I_ind_h[:, 0] / n_steps)
    plt.xlabel('t [ns]')
    plt.show()

    # Plot numerical result in 2D with particle animation
    fig = plt.figure()
    plot.get_planar_sensor_plot(fig=fig,
                                width=width,
                                pitch=pitch,
                                thickness=thickness,
                                n_pixel=n_pixel,
                                V_backplane=V_bias,
                                V_readout=V_readout,
                                potential_function=pot_descr.get_potential,
                                field_function=pot_descr.get_field,
                                depletion_function=pot_descr.get_depletion,
                                title='Planar sensor potential')
    init, animate = plot.animate_drift_diffusion(fig, pe=traj_e, ph=traj_h, dt=dt)
    ani_time = 5.
    ani = animation.FuncAnimation(fig=fig, func=animate, frames=np.arange(1, traj_h.shape[0]),
                                  interval=ani_time / traj_e.shape[0] * 1000., blit=True, init_func=init, repeat_delay=ani_time / 5.)
    plt.show()
#     plt.show()

#     x = np.arange(-2, 2)
#     y = np.arange(-2, 2)
#
#
#     fig = plt.figure()
#     ax = fig.add_subplot(111, autoscale_on=False, xlim=(-2, 2), ylim=(-2, 2))
#     ax.grid()
#
#     line, = ax.plot([], [], 'o-', lw=2)
#     time_template = 'time = %.1fs'
#     time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
#
#
#     def init():
#         line.set_data([], [])
#         time_text.set_text('')
#         return line, time_text
#
#
#     def animate(i):
#         thisx = [0, x[i]]
#         thisy = [0, y[i]]
#
#         line.set_data(thisx, thisy)
#         time_text.set_text(time_template % (i * dt))
#         return line, time_text
#
#     ani = animation.FuncAnimation(fig, animate, np.arange(1, len(y)),
#                                   interval=25, blit=True, init_func=init)

#     plt.show()
#     print p0.shape
#     print pot_descr.get_field(*p0).shape

if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    induced_current()
