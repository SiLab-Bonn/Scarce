''' Example that calculates the collected charge
    as a function of the position in a planar sensor.

    Calculates the induced charge from e-h pairs drifting
    through the silicon.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm
from scipy import integrate
from scarce import plot, solver, geometry, silicon, fields, tools, analysis, sensor

import time
from matplotlib.pyplot import twinx


def cc(bias, temperature, n_eff, t_e_trapping, t_h_trapping):
#     # Create mesh of the sensor and stores the result
#     mesh = geometry.mesh_planar_sensor(
#         n_pixel=n_pixel,
#         width=width,
#         thickness=thickness,
#         resolution=resolution)
#     # Numerically solve the Poisson equation on the mesh
#     # To create the weighting potential
#     w_potential = fields.calculate_planar_sensor_w_potential(
#         mesh=mesh,
#         width=width,
#         pitch=pitch,
#         n_pixel=n_pixel,
#         thickness=thickness)
#     min_x = float(mesh.getFaceCenters()[0, :].min())
#     max_x = float(mesh.getFaceCenters()[0, :].max())
#     nx = width * n_pixel
#     ny = thickness
#     pot_w_descr = fields.Description(w_potential,
#                                      min_x=min_x,
#                                      max_x=max_x,
#                                      min_y=0,
#                                      max_y=thickness,
#                                      nx=nx,
#                                      ny=ny,
#                                      smoothing=smoothing)
#     
#     min_x = float(mesh.getFaceCenters()[0, :].min())
#     max_x = float(mesh.getFaceCenters()[0, :].max())
# 
#     # Set um resolution grid
#     nx = width * n_pixel
#     ny = thickness
# 
#     V_bi = -silicon.get_diffusion_potential(n_eff, temperature)
#     # Numerically solve the Laplace equation on the mesh
#     potential = fields.calculate_planar_sensor_potential(
#         mesh=mesh,
#         width=width,
#         pitch=pitch,
#         n_pixel=n_pixel,
#         thickness=thickness,
#         n_eff=n_eff,
#         V_bias=bias,
#         V_readout=V_readout,
#         V_bi=V_bi)
#     pot_descr = fields.Description(potential,
#                                    min_x=min_x,
#                                    max_x=max_x,
#                                    min_y=0,
#                                    max_y=thickness,
#                                    nx=nx,
#                                    ny=ny,
#                                    smoothing=smoothing)

    # Create sensor
    pot_w_descr, pot_descr = sensor.planar_sensor(n_eff=n_eff,
                                     V_bias=bias,
                                     V_readout=V_readout,
                                     temperature=temperature,
                                     n_pixel=n_pixel,
                                     width=width,
                                     pitch=pitch,
                                     thickness=thickness,
                                     resolution=resolution,
                                     # Might have to be adjusted when changing
                                     # the geometry
                                     smoothing=smoothing
                                     )

    return analysis.get_charge_planar(width, thickness, pot_descr, pot_w_descr,
                                      t_e_trapping=t_e_trapping, t_h_trapping=t_h_trapping,
                                      grid_x=5, grid_y=5, n_pairs=20, dt=0.001,
                                      n_steps=20000, temperature=temperature), pot_descr

if __name__ == '__main__':
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")
    
    n_pixel = 9
    width = 50.
    pitch = 30.
    thickness = 250.
    V_readout = 0.

    smoothing = 0.05
    resolution = 287

    fluence = 0.
    bias = -200
    n_eff_0 = 6.475e12
    temperature=300

    n_eff = silicon.get_eff_acceptor_concentration(fluence, n_eff_0 / 1e12,
                                                         is_ntype=True,
                                                         is_oxygenated=True)[0] * 1e12
    t_e_trapping = silicon.get_trapping(fluence * 1e12, is_electron=True, paper=1)
    t_h_trapping = silicon.get_trapping(fluence * 1e12, is_electron=False, paper=1)

    (edge_x, edge_y, charge), pot_descr = cc(bias=bias,
                                            temperature=temperature,
                                            n_eff=n_eff,
                                            t_e_trapping=t_e_trapping,
                                            t_h_trapping=t_h_trapping)
    
    # Plot numerical potential result in 2D
    plot.plot_planar_sensor(width=width,
                                pitch=pitch,
                                thickness=thickness,
                                n_pixel=n_pixel,
                                V_backplane=bias,
                                V_readout=V_readout,
                                pot_func=pot_descr.get_potential,
                                field_func=pot_descr.get_field,
                                depl_func=pot_descr.get_depletion,
                                title='Planar sensor potential')

    
    # Plot collected charge map
    plt.clf()
    plt.gca().set_aspect('equal')
    plt.gca().invert_yaxis()
    cmap = cm.get_cmap('inferno')
    cmap.set_bad('white')
    cmesh = plt.pcolormesh(edge_x, edge_y, charge,
                           cmap=cmap, vmin=0, vmax=1.05)
    plt.title('Charge collection, fluence %1.2f neq_cm2' % fluence)
    plt.grid()
    cax = plt.gcf().add_axes([plt.gca().get_position().xmax, 0.1, 0.05,
                              plt.gca().get_position().ymax - plt.gca().get_position().ymin])
    plt.colorbar(cmesh, cax=cax, orientation='vertical')
    plt.grid()
    plt.savefig('CC_%d_%d.pdf' % (fluence, bias), layout='tight')
    plt.show()
