r"""The solver module provides numerical ODE solver functions.

These functions are either implemented here or provide an interface
to other solvers.
"""

import fipy
import numpy as np
import matplotlib.pyplot as plt
import progressbar

from scarce import silicon


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
    tolerance = 1e-15
    iterations = 10000

    equation.solve(var=var,
                   solver=solver(tolerance=tolerance, iterations=iterations),
                   **kwargs)


class DriftDiffusionSolver(object):

    ''' Solve the drift-diffusion equation for pseudo particles.

        Solving via euler forward difference.
    '''

    def __init__(self, pot_descr, pot_w_descr,
                 T=300):
        '''
        Parameters
        ----------
        pot_descr : fields.Description
            Description object to describe the potential

        pot_w_descr : fields.Description
            Description object to describe the weighting potential

        T : number
            Temperatur in Kelvin

        Notes
        -----

        '''

        self.pot_descr = pot_descr
        self.pot_w_descr = pot_w_descr

        self.T = T

    def solve(self, p0, q0, dt, n_steps=4000):
        ''' Solve the drift diffusion equation for quasi partciles and calculates
            the total induced current.

            Parameters
            ----------
            p0 : array_like, shape = (n, 2)
                n positions of e-h-pairs in x, y at t = 0

            q0 : array_like, shape = (n, )
                n charges of quasi particles at t = 0

            dt : number
                Time step in ns. Influences presicion.

            n_steps : number
                Number of steps in time. Has to be large enough that all
                particles reach the readout electrode.
        '''

        # Start positions
        p_e, p_h = p0, p0.copy()

        # Result arrays initialized to NaN
        traj_e = np.empty(shape=(n_steps, p_e.shape[0], p_e.shape[1]))
        traj_h = np.empty(shape=(n_steps, p_h.shape[0], p_h.shape[1]))
        I_ind_e = np.empty(shape=(n_steps, p_e.shape[1]))
        I_ind_h = np.empty(shape=(n_steps, p_h.shape[1]))
        traj_e[:], traj_h[:] = np.nan, np.nan
        I_ind_e[:], I_ind_h[:] = np.nan, np.nan

        progress_bar = progressbar.ProgressBar(
            widgets=['', progressbar.Percentage(), ' ',
                     progressbar.Bar(marker='*', left='|', right='|'), ' ',
                     progressbar.AdaptiveETA()],
            maxval=n_steps, term_width=80)

        progress_bar.start()

        def in_boundary(description, x, y):
            ''' Checks if the particles are still in the bulk and should be propagated
            '''

            sel_x = np.logical_and(x >= description.min_x,
                                   x <= description.max_x)
            sel_y = np.logical_and(y >= description.min_y,
                                   y <= description.max_y)
            return np.logical_and(sel_x, sel_y)

        sel_e = np.ones(p_e.shape[1], dtype=np.bool)
        sel_h = np.ones(p_h.shape[1], dtype=np.bool)

        for step in range(n_steps):
            # Check if all particles out of boundary
            if not np.any(sel_e) and not np.any(sel_h):
                break

            # Electric field in V/cm
            E_e = self.pot_descr.get_field(p_e[0, sel_e], p_e[1, sel_e]) * 1e4
            E_h = self.pot_descr.get_field(p_h[0, sel_h], p_h[1, sel_h]) * 1e4

            # Mobility in cm2 / Vs
            mu_e = silicon.get_mobility(np.sqrt(E_e[0] ** 2 + E_e[1] ** 2),
                                        temperature=self.T, is_electron=True)
            mu_h = silicon.get_mobility(np.sqrt(E_h[0] ** 2 + E_h[1] ** 2),
                                        temperature=self.T, is_electron=False)

            # Velocity in cm / s
            v_e, v_h = - E_e * mu_e, E_h * mu_h

            # Calculate induced current
            if np.any(p_e[0, sel_e]):  # Only if electrons are still drifting
                # Weighting field in V/um
                W_e = self.pot_w_descr.get_field(p_e[0, sel_e],
                                                 p_e[1, sel_e])
                # Induced current in C/s
                I_ind_e[step, sel_e] = (W_e[0, 0] * v_e[0, 0] +
                                        W_e[1, 0] * v_e[1, 0]) \
                    - q0[sel_e] * dt * 1e-5

            if np.any(p_h[0, sel_h]):  # Only if holes are still drifting
                # Weighting field in V/um
                W_h = self.pot_w_descr.get_field(p_h[0, sel_h],
                                                 p_h[1, sel_h])
                # Induced current in C/s
                I_ind_h[step, sel_h] = (W_h[0, 0] * v_h[0, 0] +
                                        W_h[1, 0] * v_h[1, 0]) * \
                    q0[sel_h] * dt * 1e-5

            # Position change in um
            d_p_e, d_p_h = v_e * dt * 1e-5, v_h * dt * 1e-5

            # Update position
            p_e[:, sel_e] = p_e[:, sel_e] + d_p_e
            p_h[:, sel_h] = p_h[:, sel_h] + d_p_h

            # Check boundaries and update selection
            sel_e = in_boundary(description=self.pot_descr,
                                x=p_e[0, :], y=p_e[1, :])
            sel_h = in_boundary(description=self.pot_descr,
                                x=p_h[0, :], y=p_h[1, :])
            p_e[:, ~sel_e] = np.nan
            p_h[:, ~sel_h] = np.nan

            traj_e[step] = p_e
            traj_h[step] = p_h

            progress_bar.update(step)
        progress_bar.finish()

        return traj_e, traj_h, I_ind_e, I_ind_h
