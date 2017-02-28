r"""Helper functions for complex analysis """

import numpy as np
from scipy import integrate

from scarce import solver


def get_charge_planar(width, thickness, pot_descr, pot_w_descr, t_e_trapping=0., t_h_trapping=0., grid_x=5, grid_y=5, n_pairs=10, dt=0.001, n_steps=25000, temperature=300):
        ''' Calculate the collected charge in one planar pixel

            Charge is given as a 2d map depending on the start postitions of the e-h pairs.

            Parameters
            ----------
            width: number
                Pixel width in um
            thickness: number
                Pixel thickness in um
            pot_descr, pot_w_descr: scarce.fields.Description
                Solution for the drift/weightning potential.
            grid_x, grid_y: number
                Grid spacing in um
            n_pairs: number
                of pseudo e-h pairs per grid point
            dt: float
                Time step in simulation in ns. Should be 1 ps to give reasonable diffusion
            n_steps: int
                Time steps to simulate
        '''

        # Number of x/y bins
        x_bins = int(width / grid_x)
        y_bins = int(thickness / grid_y)
        # Bin positions
        range_x = (-width / 2., width / 2.)
        range_y = (0, thickness)
        # Create e-h pairs in the pixel, avoid charge carriers on boundaries
        # e.g. x = -width / 2 or y = 0
        xx, yy = np.meshgrid(np.linspace(range_x[0] + grid_x / 2.,
                                         range_x[1] - grid_x / 2., x_bins),
                             np.repeat(np.linspace(range_y[0] + grid_y / 2.,
                                                   range_y[1] - grid_y / 2., y_bins),
                                       n_pairs),  # 10 e-h per position
                             sparse=False)  # All combinations of x / y
        # Start positions
        p0 = np.array([xx.ravel(), yy.ravel()])  # Position [um]

        # Initial charge set to 1
        q_start = 1.
        q0 = np.ones(p0.shape[1]) * q_start
        # Needed for histograming, numerical accuracy demands > 1
        q_max = q_start * 1.05

        t = np.linspace(0, n_steps * dt, 1000)

        dd = solver.DriftDiffusionSolver(pot_descr, pot_w_descr,
                                         T=temperature, diffusion=True,
                                         t_e_trapping=t_e_trapping, t_h_trapping=t_h_trapping, save_frac=50)
        traj_e, traj_h, I_ind_e, I_ind_h, T, _, Q_ind_e_tot, Q_ind_h_tot = dd.solve(p0, q0, dt, n_steps,
                                                                                    multicore=True)

        # Trajectory at t=0 is start position
        pos_0 = traj_e[0]

        # Interpolate data to fixed time points for easier plotting
    #     I_ind_e = tools.time_data_interpolate(T, I_ind_e, t, axis=0, fill_value=0.)
    #     I_ind_h = tools.time_data_interpolate(T, I_ind_h, t, axis=0, fill_value=0.)
        I_ind_e[np.isnan(I_ind_e)] = 0.
        I_ind_h[np.isnan(I_ind_h)] = 0.
        Q_ind_e = integrate.cumtrapz(I_ind_e, T, axis=0, initial=0)
        Q_ind_h = integrate.cumtrapz(I_ind_h, T, axis=0, initial=0)

        # Last index with data (time != nan)
    #     index = np.nanargmax(T, axis=0)
    #     y = np.indices(index.shape)
        # Last recorded integrated charge is total induced charge
    #     q_ind = Q_ind_e[index, y][0] + Q_ind_h[index, y][0]
        q_ind = Q_ind_e_tot + Q_ind_h_tot

        # Histogram charge per start position
        data = np.vstack((pos_0[0], pos_0[1], q_ind)).T
        n_bins_c = 200  # Number of charge bins
        H, edges = np.histogramdd(sample=data,
                                  bins=(x_bins, y_bins, n_bins_c),
                                  range=((range_x[0], range_x[1]),
                                         (range_y[0], range_y[1]),
                                         (0., q_max)))
        # Result hist
        charge_pos = np.zeros(shape=(x_bins, y_bins))
        sel = (np.sum(H, axis=2) != 0)
        weights = (edges[2][:-1] + edges[2][1:]) / 2.
        charge_pos[sel] = np.average(H, axis=2,
                                     weights=weights)[sel] * weights.sum() / np.sum(H, axis=2)[sel]
        edges_x = (edges[0][:-1] + edges[0][1:]) / 2.
        edges_y = (edges[1][:-1] + edges[1][1:]) / 2.

#         for xi, yi in zip(*np.where(np.logical_and(charge_pos > 0.1,
#                                                    charge_pos < 0.9))
#                           ):
#             print edges_x[xi], edges_y[yi], charge_pos[xi, yi]
#             plt.clf()
#             plt.bar(weights, H[xi, yi], width=np.diff(weights)[0])
#             plt.show()
#             plt.clf()
#             sel = np.logical_and(pos_0[0] == edges_x[xi],
#                                  pos_0[1] == edges_y[yi])
#             plt.plot(T[:, sel], Q_ind_e[:, sel] + Q_ind_h[:, sel])
#             for c in weights[H[xi, yi].astype(np.bool)]:
#                 plt.plot(plt.xlim(), [c, c])
#             plt.show()
#             plt.clf()
#             plt.plot(
#                 T[:, sel], traj_e[:, 0, sel], '-.', linewidth=1, label='e_x')
#             plt.plot(
#                 T[:, sel], traj_e[:, 1, sel], '--', linewidth=1, label='e_y')
#             plt.plot(
#                 T[:, sel], traj_h[:, 0, sel], '-.', linewidth=1, label='h_x')
#             plt.plot(
#                 T[:, sel], traj_h[:, 1, sel], '--', linewidth=1, label='h_y')
#             plt.legend(loc=2)
#             plt.show()
#             break

        return edges[0], edges[1], charge_pos.T
