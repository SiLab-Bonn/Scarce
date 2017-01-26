r"""The solver module provides numerical ODE solver functions.

These functions are either implemented here or provide an interface
to other solvers.
"""

import fipy
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count
from functools import partial
import progressbar
import logging
from scarce import silicon, tools


def solve(var, equation, **kwargs):
    ''' Interface to the fipy solver used for the 2d poisson equation.

    Parameters
    ----------
    var : fipy.CellVariable
          Fipy meshed cell variables to solve the equation for

    equation : fipy equation
       e.g. fipy.DiffusionTerm() = 0

    kwargs : kwargs
        Arguments of fipy equation.solve(kwargs)

    Notes
    -----
    A linear LU solver is forced here, since otherwise pysparse based solver
    do not converge properly.
    '''

    # Force a linear solver
    # The others do not always converge (?!)
    solver = fipy.solvers.LinearLUSolver

    # Reasonable and more strict convergence criteria
    tolerance = 1e-20
    iterations = 10000

    equation.solve(var=var,
                   solver=solver(tolerance=tolerance, iterations=iterations),
                   **kwargs)


class DriftDiffusionSolver(object):

    ''' Solve the drift-diffusion equation for pseudo particles.

        Solving via euler forward difference.
    '''

    def __init__(self, pot_descr, pot_w_descr,
                 T=300, geom_descr=None, diffusion=True,
                 t_e_trapping=0., t_h_trapping=0.):
        '''
        Parameters
        ----------
        pot_descr : fields.Description
            Description object to describe the potential

        pot_w_descr : fields.Description
            Description object to describe the weighting potential

        T : number
            Temperatur in Kelvin

        geom_descr : scarce.geometry.SensorDescription3D
            Describes the 3D sensor geometry

        t_e_trapping : number
            Trapping time for electrons in ns

        t_h_trapping : number
            Trapping time for holes in ns

        Notes
        -----

        '''

        self.pot_descr = pot_descr
        self.pot_w_descr = pot_w_descr

        self.T = T  # Temperatur in Kelvin
        self.t_e_trapping = t_e_trapping  # Trapping time in ns
        self.t_h_trapping = t_h_trapping  # Trapping time in ns
        self.geom_descr = geom_descr
        self.diffusion = diffusion

    def solve(self, p0, q0, dt, n_steps=4000):
        ''' Solve the drift diffusion equation for quasi partciles and calculates
            the total induced current.

            Parameters
            ----------
            p0 : array_like, shape = (n, 2)
                n positions of e-h-pairs in x, y at t = 0

            q0 : array_like, shape = (n, )
                Charge of quasi particles at t = 0.
                To give in electrons is a good choise.

            dt : number
                Time step in ns. Influences presicion. Should be <= 0.001 (ps).

            n_steps : number
                Number of steps in time. Has to be large enough that all
                particles reach the readout electrode.
        '''

        if dt > 0.001:
            logging.warning('A time step > 1 ps result in wrong diffusion')

        # Split data and into cores - 1 slices
        pool = Pool()
        n_slices = cpu_count() - 1
        slices_p0 = np.array_split(p0, n_slices, axis=1)
        slices_q0 = np.array_split(q0, n_slices, axis=0)

        logging.info('Calculate drift diffusion on %d cores', n_slices)

        jobs = []
        for index in range(n_slices):
            job = tools.apply_async(pool=pool,
                                    fun=_solve_dd,
                                    p0=slices_p0[index],
                                    q0=slices_q0[index],
                                    n_steps=n_steps,
                                    dt=dt,
                                    geom_descr=self.geom_descr,
                                    pot_w_descr=self.pot_w_descr,
                                    pot_descr=self.pot_descr,
                                    T=self.T,
                                    diffusion=self.diffusion,
                                    t_e_trapping=self.t_e_trapping,
                                    t_h_trapping=self.t_h_trapping
                                    )
            jobs.append(job)

        # Gather results
        results = []
        for job in jobs:
            results.append(job.get())

        # Merge results
        traj_e = np.concatenate([i[0] for i in results], axis=2)
        traj_h = np.concatenate([i[1] for i in results], axis=2)
        I_ind_e = np.concatenate([i[2] for i in results], axis=1)
        I_ind_h = np.concatenate([i[3] for i in results], axis=1)

        pool.close()
        pool.join()

        return traj_e, traj_h, I_ind_e, I_ind_h


def test(**kwarg):
    pass

# Drift diffussion iteration loop helper functions


def _in_boundary(geom_descr, pot_descr, x, y):
    ''' Checks if the particles are still in the bulk and should be propagated
    '''

    sel_x = np.logical_and(x >= pot_descr.min_x,
                           x <= pot_descr.max_x)
    sel_y = np.logical_and(y >= pot_descr.min_y,
                           y <= pot_descr.max_y)

    sel = np.logical_and(sel_x, sel_y)

    if geom_descr:
        sel_col = geom_descr.position_in_column(x, y, incl_sides=True)
        sel = np.logical_and(sel, ~sel_col)

    return sel


def _correct_boundary(geom_descr, pot_descr, x, y, is_electron):
    ''' Checks if the particles are still in the bulk and should be propagated
    '''

    if not geom_descr:  # planar sensor
        if is_electron:
            x[x >= pot_descr.max_x] = pot_descr.max_x
        else:
            x[x <= 0.] = 0.


def _solve_dd(p0, q0, n_steps, dt, geom_descr, pot_w_descr, pot_descr, T,
              diffusion, t_e_trapping, t_h_trapping):
    # E-h pairs Start positions
    p_e, p_h = p0.copy(), p0.copy()

    # Result arrays initialized to NaN
    traj_e = np.empty(shape=(n_steps, p_e.shape[0], p_e.shape[1]))
    traj_h = np.empty(shape=(n_steps, p_h.shape[0], p_h.shape[1]))
    I_ind_e = np.zeros(shape=(n_steps, p_e.shape[1]))
    I_ind_h = np.zeros(shape=(n_steps, p_h.shape[1]))
    traj_e[:], traj_h[:] = np.nan, np.nan

#     progress_bar = progressbar.ProgressBar(
#         widgets=['', progressbar.Percentage(), ' ',
#                  progressbar.Bar(marker='*', left='|', right='|'), ' ',
#                  progressbar.AdaptiveETA()],
#         maxval=n_steps, term_width=80)
#
#     progress_bar.start()

    sel_e = _in_boundary(geom_descr, pot_descr=pot_descr,
                         x=p_e[0, :], y=p_e[1, :])
    sel_h = _in_boundary(geom_descr, pot_descr=pot_descr,
                         x=p_h[0, :], y=p_h[1, :])

    for step in range(n_steps):
        # Store position in trajectory arrays
        traj_e[step] = p_e
        traj_h[step] = p_h
        # Check if all particles out of boundary
        if not np.any(sel_e) and not np.any(sel_h):
            break  # Stop loop to safe time

        # Electric field in V/cm
        E_e = pot_descr.get_field(
            p_e[0, sel_e], p_e[1, sel_e]) * 1e4
        E_h = pot_descr.get_field(
            p_h[0, sel_h], p_h[1, sel_h]) * 1e4

        # Mobility in cm2 / Vs
        mu_e = silicon.get_mobility(np.sqrt(E_e[0] ** 2 + E_e[1] ** 2),
                                    temperature=T, is_electron=True)
        mu_h = silicon.get_mobility(np.sqrt(E_h[0] ** 2 + E_h[1] ** 2),
                                    temperature=T, is_electron=False)

        # Drift velocity in cm / s
        v_e, v_h = - E_e * mu_e, E_h * mu_h

        # Add diffusion velocity
        if diffusion:
            # Calculate absolute thermal velocity
            v_th_e = silicon.get_thermal_velocity(temperature=T,
                                                  is_electron=True)
            v_th_h = silicon.get_thermal_velocity(temperature=T,
                                                  is_electron=False)
            # Create thermal velocity distribution
            # From: IEEE VOL. 56, NO. 3, JUNE 2009
            v_th_e *= np.log(np.abs(
                1. / (1. - np.random.uniform(size=v_e.shape[1]))))
            v_th_h *= np.log(np.abs(
                1. / (1. - np.random.uniform(size=v_h.shape[1]))))
            # Calculate random direction in x, y
            # Uniform random number 0 .. 2 Pi
            eta = np.random.uniform(0., 2. * np.pi, size=v_e.shape[1])
            direction_e = np.array([np.cos(eta), np.sin(eta)])
            eta = np.random.uniform(0., 2. * np.pi, size=v_h.shape[1])
            direction_h = np.array([np.cos(eta), np.sin(eta)])

            v_th_e = v_th_e[np.newaxis, :] * direction_e
            v_th_h = v_th_h[np.newaxis, :] * direction_h

            v_e += v_th_e
            v_h += v_th_h

        # Calculate induced current
        # Only if electrons are still drifting
        if np.any(sel_e):
            # Weighting field in V/um
            W_e = pot_w_descr.get_field(p_e[0, sel_e], p_e[1, sel_e])

            # Induced charge in C/s, Q = E_w * v * q * dt
            dQ_e = (W_e[0] * v_e[0] + W_e[1] * v_e[1]) * \
                - q0[sel_e] * dt * 1e-5

            # Reduce induced charge due to trapping
            if t_e_trapping:
                dQ_e *= np.exp(-dt * step / t_e_trapping)

            # Induced current
            I_ind_e[step, sel_e] = dQ_e / dt

        if np.any(sel_h):  # Only if holes are still drifting
            # Weighting field in V/um
            W_h = pot_w_descr.get_field(p_h[0, sel_h], p_h[1, sel_h])

            # Induced charge in C/s, Q = E_w * v * q * dt
            dQ_h = (W_h[0] * v_h[0] + W_h[1] * v_h[1]) * \
                q0[sel_h] * dt * 1e-5

            # Reduce induced charge due to trapping
            if t_h_trapping:
                dQ_h *= np.exp(-dt * step / t_e_trapping)

            # Induced current
            I_ind_h[step, sel_h] = dQ_h / dt

        # Position change in um
        d_p_e, d_p_h = v_e * dt * 1e-5, v_h * dt * 1e-5

        # Update position
        p_e[:, sel_e] = p_e[:, sel_e] + d_p_e
        p_h[:, sel_h] = p_h[:, sel_h] + d_p_h

        # Correct boundaries (e.g. leaving sensor due to diffusion)
        _correct_boundary(geom_descr, pot_descr,
                          x=p_e[0, :], y=p_e[1, :], is_electron=True)
        _correct_boundary(geom_descr, pot_descr,
                          x=p_h[0, :], y=p_h[1, :], is_electron=False)

        # Check boundaries and update selection
        sel_e = _in_boundary(geom_descr, pot_descr=pot_descr,
                             x=p_e[0, :], y=p_e[1, :])
        sel_h = _in_boundary(geom_descr, pot_descr=pot_descr,
                             x=p_h[0, :], y=p_h[1, :])
        p_e[:, ~sel_e] = np.nan
        p_h[:, ~sel_h] = np.nan

#         progress_bar.update(step)
#     progress_bar.finish()

    return traj_e, traj_h, I_ind_e, I_ind_h
