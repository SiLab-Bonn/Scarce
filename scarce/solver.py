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
                 t_e_trapping=0., t_h_trapping=0., 
                 t_e_t1=0., t_h_t1=0.,
                 t_r=0., save_frac=20):
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

        t_e_t1 : number
            E-Field dependent trapping time in ns / (V / um)

        t_h_t1 : number
            E-Field dependent trapping time in ns / (V / um)

        t_r : number
            Peaking time of CSA

        save_frac : number
            Fraction of time steps to save for each e-h pair.
            E.g.: 100 means that the position and current is saved for every
            100th time step. Be aware that the sampling might be too low to get
            the maxima correctly. If you want to be sure to get the correct
            single e-h pair values save_frac should be set to 1. Usually one
            is interessted in the total current. This value is independent of
            save_frac but not available for each e-h pair.

        Notes
        -----

        '''

        self.pot_descr = pot_descr
        self.pot_w_descr = pot_w_descr

        self.T = T  # Temperatur in Kelvin
        self.t_e_trapping = t_e_trapping  # Trapping time in ns
        self.t_h_trapping = t_h_trapping  # Trapping time in ns
        self.t_e_t1 = t_e_t1  # Trapping time in ns
        self.t_h_t1 = t_h_t1  # Trapping time in ns
        self.t_r = t_r  # Rise time of Amplifier
        self.geom_descr = geom_descr
        self.diffusion = diffusion

        self.save_frac = save_frac

    def solve(self, p0, q0, dt, n_steps, multicore=True):
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

        if not multicore:
            # E-h pairs Start positions
            p_e_0, p_h_0 = p0.copy(), p0.copy()
            return _solve_dd(p_e_0,
                             p_h_0,
                             q0=q0,
                             n_steps=n_steps,
                             dt=dt,
                             geom_descr=self.geom_descr,
                             pot_w_descr=self.pot_w_descr,
                             pot_descr=self.pot_descr,
                             temp=self.T,
                             diffusion=self.diffusion,
                             t_e_trapping=self.t_e_trapping,
                             t_h_trapping=self.t_h_trapping,
                             t_e_t1=self.t_e_t1,
                             t_h_t1=self.t_h_t1,
                             t_r=self.t_r,
                             save_frac=self.save_frac)
            
        n_slices = min(4, cpu_count() - 1)
        pool = Pool(n_slices)

        # Split data and into cores - 1 slices, max 4
        slices_p_e_0 = np.array_split(p0, n_slices, axis=1)
        slices_p_h_0 = np.array_split(p0, n_slices, axis=1)
        slices_q0 = np.array_split(q0, n_slices, axis=0)

        logging.info('Calculate drift diffusion on %d cores', n_slices)

        jobs = []
        for index in range(n_slices):
            job = tools.apply_async(pool=pool,
                                    fun=_solve_dd,
                                    p_e_0=slices_p_e_0[index],
                                    p_h_0=slices_p_h_0[index],
                                    q0=slices_q0[index],
                                    n_steps=n_steps,
                                    dt=dt,
                                    geom_descr=self.geom_descr,
                                    pot_w_descr=self.pot_w_descr,
                                    pot_descr=self.pot_descr,
                                    temp=self.T,
                                    diffusion=self.diffusion,
                                    t_e_trapping=self.t_e_trapping,
                                    t_h_trapping=self.t_h_trapping,
                                    t_e_t1=self.t_e_t1,
                                    t_h_t1=self.t_h_t1,
                                    t_r=self.t_r,
                                    save_frac=self.save_frac
                                    )
            jobs.append(job)

        # Gather results
        results = []
        for job in jobs:
            results.append(job.get())

        pool.close()
        pool.join()
        
        del pool

        # Merge results
        I_ind_tot = np.zeros(shape=(n_steps,))
        traj_e = np.concatenate([i[0] for i in results], axis=2)
        traj_h = np.concatenate([i[1] for i in results], axis=2)
        I_ind_e = np.concatenate([i[2] for i in results], axis=1)
        I_ind_h = np.concatenate([i[3] for i in results], axis=1)
        T = np.concatenate([i[4] for i in results], axis=1)
        for i in results:
            I_ind_tot += i[5]
        Q_ind_tot_e = np.concatenate([i[6] for i in results], axis=0)
        Q_ind_tot_h = np.concatenate([i[7] for i in results], axis=0)

        return traj_e, traj_h, I_ind_e, I_ind_h, T, I_ind_tot, Q_ind_tot_e, Q_ind_tot_h


# Drift diffusion iteration loop helper functions
def _in_boundary(geom_descr, pot_descr, x, y, sel):
    ''' Checks if the particles are still in the bulk and should be propagated
    '''

    sel_x = np.logical_and(x[sel] >= pot_descr.min_x,
                           x[sel] <= pot_descr.max_x)
    sel_y = np.logical_and(y[sel] >= pot_descr.min_y,
                           y[sel] <= pot_descr.max_y)

    old_sel = sel.copy()

    sel[sel] = np.logical_and(sel_x, sel_y)

    if geom_descr:
        sel_col = geom_descr.position_in_column(x, y, incl_sides=True)
        sel = np.logical_and(sel, ~sel_col)

    return sel, np.where(old_sel != sel)[0]


def _correct_boundary(geom_descr, pot_descr, x, y, sel, is_electron):
    ''' Checks if the particles are still in the bulk and should be propagated
    '''

    if not geom_descr:  # planar sensor
        if is_electron:
            y[sel][y[sel] > pot_descr.max_y] = pot_descr.max_y
        else:
            y[sel][y[sel] < 0.] = 0.


def _solve_dd(p_e_0, p_h_0, q0, n_steps, dt, geom_descr, pot_w_descr,
              pot_descr, temp, diffusion, t_e_trapping, t_h_trapping,
              t_e_t1, t_h_t1, t_r, save_frac):
    p_e, p_h = p_e_0, p_h_0

    # Result arrays initialized to NaN
    n_store = int(n_steps / save_frac)  # Steps to store

    max_step_size = n_steps / n_store * 10
    # Different store time step for each e-h pair
    T = np.full(shape=(n_store, p_e.shape[1]),
                fill_value=np.nan, dtype=np.float32)
    # Stored trajectory for each eh pair
    traj_e = np.full(shape=(n_store, p_e.shape[0], p_e.shape[1]),
                     fill_value=np.nan)
    traj_h = np.full(shape=(n_store, p_h.shape[0], p_h.shape[1]),
                     fill_value=np.nan)
    # Stored induced charge for each eh pair
    I_ind_e = np.zeros(shape=(n_store, p_e.shape[1]))
    I_ind_h = np.zeros_like(I_ind_e)
    # Helper array(s) of actual step index and next step index
    i_step = np.zeros(p_e.shape[1], dtype=np.int)  # Result array indeces
    next_step = np.zeros_like(i_step)  # Next time step to store
    # Summed induced charge/current with every time step
    I_ind_tot = np.zeros(shape=(n_steps))
    # Total induced charge
    Q_ind_tot_e = np.zeros(shape=(p_e.shape[1]))
    Q_ind_tot_h = np.zeros(shape=(p_h.shape[1]))

    def add_diffusion(v_e, v_h):
        # Calculate absolute thermal velocity
        v_th_e = silicon.get_thermal_velocity(temperature=temp,
                                              is_electron=True)
        v_th_h = silicon.get_thermal_velocity(temperature=temp,
                                              is_electron=False)
        # Create thermal velocity distribution
        # From: The Atomistic Simulation of Thermal Diffusion
        # and Coulomb Drift in Semiconductor Detectors
        # IEEE VOL. 56, NO. 3, JUNE 2009
        v_th_e *= np.sqrt(2. / 3. * np.log(np.abs(
            1. / (1. - np.random.uniform(size=v_e.shape[1])))))
        v_th_h *= np.sqrt(2. / 3. * np.log(np.abs(
            1. / (1. - np.random.uniform(size=v_h.shape[1])))))
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

        return v_e, v_h

    def cal_step_size(t, Q_ind_tot, dt, dydt, q_max, i_step, n_store):
        ''' Calculates the step size from the actual data

        It is assumed that the actual slope stays constant.
        The remaining time distance is calculated and the step size is
        adjusted to fit the remaining steps.
        '''

        step_size = np.zeros_like(q_max)

        # All steps are used, mark as done
        sel_done = n_store - i_step - 1 <= 0

        # All storage spaces used up
        if step_size[~sel_done].size == 0:
            return step_size.astype(np.int)

        # Set next step to NaN if storing is done
        step_size[sel_done] = np.NaN

        # Calculate remaining x distance to cover
        try:
            # Case: increasing function (assume value = q_max at tmax)
            t_max_exp = ((q_max - Q_ind_tot) / dydt + t)[~sel_done]
            # Case: decreasing function (assume value = 0 at tmax)
            t_max_exp[dydt[~sel_done] < 0] = ((-q_max - Q_ind_tot) /
                                              dydt + t)[~sel_done][dydt[~sel_done] < 0]
            # Correct expected time larger than simulation time
            t_max_exp[t_max_exp > n_steps * dt] = n_steps * dt
            # Remaining time to cover
            t_left = t_max_exp - t
        except (IndexError, ValueError):
            logging.error('q_max.shape %s', str(q_max.shape))
            logging.error('sel_done.shape %s', str(sel_done.shape))
            logging.error('Q_ind_tot.shape %s', str(Q_ind_tot.shape))
            logging.error('dydt.shape %s', str(dydt.shape))
            logging.error('t_max_exp.shape %s', str(t_max_exp.shape))
            logging.error('dydt[~sel_done].shape %s', str(dydt[~sel_done].shape))
            logging.error('(dydt[~sel_done] < 0).shape %s', str((dydt[~sel_done] < 0).shape))
            logging.error('(-q_max - Q_ind_tot) / dydt + t).shape %s', str(((-q_max - Q_ind_tot) / dydt + t).shape))
            logging.error('(-q_max - Q_ind_tot) / dydt + t)[~sel_done].shape %s', str(((-q_max - Q_ind_tot) / dydt + t)[~sel_done].shape))

            raise

        # Calculate the step size
        step_size[~sel_done] = t_left / dt / (n_store - i_step[~sel_done] - 1)

        # Limit step size to max_step_size
        # Needed for slope direction changing functions
        sel_lim_max = step_size[~sel_done] > max_step_size
        step_size[~sel_done][sel_lim_max] = max_step_size
        # Minimum step size = 1
        step_size[np.logical_and(~sel_done, step_size < 1.)] = 1

        step_size = step_size.astype(np.int)

        return step_size

    def store_if_needed(step, next_step, Q_ind_tot_e, Q_ind_tot_h, dt, dQ_e,
                        dQ_h, T, I_ind_e, I_ind_h, p_e, p_h, q_max,
                        i_step, n_store, sel_e, sel_h):
        ''' Checks if charge carriers value needs to be stored and returns
                next storage time step
        '''

        t = dt * step  # Actual time step

        # Select charges that need storing for full array (1 entry per e-h)
        store = step == next_step

        # Select charges that need storing and are not fully propagated yet
        # for full array (1 entry per e-h)
        store_e = np.logical_and(store, sel_e)
        store_h = np.logical_and(store, sel_h)

        # Select charges that need storing for reduced array (e-h that need
        # propagating
        s_e = store_e[sel_e]
        s_h = store_h[sel_h]

        # All carrier stored
        if not np.any(store):
            return

        # All carriers needing storing are fully drifted
        if not np.any(s_e) and not np.any(s_h):
            return

        # Set data of actual time step
        T[i_step[store], store] = t

        # Calculate induced charge for integrated time steps T
        DT_e = T[i_step[store_e], store_e] - T[i_step[store_e] - 1, store_e]
        DT_h = T[i_step[store_h], store_h] - T[i_step[store_h] - 1, store_h]
        DT_e[np.isnan(DT_e)] = dt  # First storing
        DT_h[np.isnan(DT_h)] = dt  # First storing
        I_ind_e[i_step[store_e], store_e] = dQ_e_step[store_e] / DT_e
        I_ind_h[i_step[store_h], store_h] = dQ_h_step[store_h] / DT_h

        # Store data
#         I_ind_e[i_step[store_e], store_e] = dQ_e[s_e] / dt
#         I_ind_h[i_step[store_h], store_h] = dQ_h[s_h] / dt

        traj_e[i_step[store_e], :, store_e] = p_e[:, store_e].T
        traj_h[i_step[store_h], :, store_h] = p_h[:, store_h].T

#         if np.max(i_step) == 145:
#             plt.plot(T[:, 0], I_ind_e[:, 0], '.')
#             plt.plot(T[:, 0], I_ind_h[:, 0], '.')
#             plt.show()

        # Calculate step size as the minimum of the e and h step size
        d_step = np.zeros_like(next_step)
        d_step_e = np.zeros_like(d_step)
        d_step_h = np.zeros_like(d_step)
        if np.any(store_e):
            d_step_e[store_e] = cal_step_size(t, Q_ind_tot=Q_ind_tot_e[store_e], dt=dt,
                                              dydt=I_ind_e[i_step[store_e], store_e],
                                              q_max=q_max[store_e],
                                              i_step=i_step[store_e],
                                              n_store=n_store)
            d_step[store_e] = d_step_e[store_e]

        if np.any(store_h):
            d_step_h[store_h] = cal_step_size(t, Q_ind_tot=Q_ind_tot_h[store_h], dt=dt,
                                              dydt=I_ind_h[i_step[store_h], store_h],
                                              q_max=q_max[store_h],
                                              i_step=i_step[store_h],
                                              n_store=n_store)
            d_step[store_h] = d_step_h[store_h]
        sel = np.logical_and(store_e, store_h)
        d_step[sel] = np.minimum(d_step_e[sel], d_step_h[sel])

        next_step += d_step

        # Increase storage hists indeces
        i_step[store] += 1

        i_step[i_step >= T.shape[0]] = T.shape[0] - 1

        # Reset tmp. variable
        dQ_e_step[store_e] = 0.
        dQ_h_step[store_h] = 0.

    # Tmp. variables to store the total induced charge per save step
    # Otherwise the resolution of induced current calculation
    # is reduced
    dQ_e_step = np.zeros(shape=p_e.shape[1])
    dQ_h_step = np.zeros(shape=p_h.shape[1])

    progress_bar = progressbar.ProgressBar(
        widgets=['', progressbar.Percentage(), ' ',
                 progressbar.Bar(marker='*', left='|', right='|'), ' ',
                 progressbar.AdaptiveETA()],
        maxval=n_steps, term_width=80)

    progress_bar.start()

    sel_e, _ = _in_boundary(geom_descr, pot_descr=pot_descr,
                            x=p_e[0, :], y=p_e[1, :],
                            sel=np.ones(p_e.shape[1], dtype=np.bool))
    sel_h, _ = _in_boundary(geom_descr, pot_descr=pot_descr,
                            x=p_h[0, :], y=p_h[1, :],
                            sel=np.ones(p_h.shape[1], dtype=np.bool))

    for step in range(n_steps):
        # Check if all particles out of boundary
        if not np.any(sel_h) and not np.any(sel_e):
            break  # Stop loop to safe time

        # Electric field in V/cm
        E_e = pot_descr.get_field(
            p_e[0, sel_e], p_e[1, sel_e]) * 1e4
        E_h = pot_descr.get_field(
            p_h[0, sel_h], p_h[1, sel_h]) * 1e4

        # Mobility in cm2 / Vs
        mu_e = silicon.get_mobility(np.sqrt(E_e[0] ** 2 + E_e[1] ** 2),
                                    temperature=temp, is_electron=True)
        mu_h = silicon.get_mobility(np.sqrt(E_h[0] ** 2 + E_h[1] ** 2),
                                    temperature=temp, is_electron=False)

        # Drift velocity in cm / s
        v_e, v_h = - E_e * mu_e, E_h * mu_h

        if diffusion:
            v_e, v_h = add_diffusion(v_e, v_h)

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
                t_e = t_e_trapping + t_e_t1 * np.sqrt(E_e[0]**2 + E_e[1]**2) * 1e-4
                dQ_e *= np.exp(-dt * step / t_e)

            if t_r:
                dQ_e *= np.exp(-dt * step / (t_r / 2.2))

            dQ_e_step[sel_e] += dQ_e

            Q_ind_tot_e[sel_e] += dQ_e

            I_ind_tot[step] += dQ_e.sum() / dt

        if np.any(sel_h):  # Only if holes are still drifting
            # Weighting field in V/um
            W_h = pot_w_descr.get_field(p_h[0, sel_h], p_h[1, sel_h])

            # Induced charge in C/s, Q = E_w * v * q * dt
            dQ_h = (W_h[0] * v_h[0] + W_h[1] * v_h[1]) * \
                q0[sel_h] * dt * 1e-5

            # Reduce induced charge due to trapping
            if t_h_trapping:
                t_h = t_h_trapping + t_h_t1 * np.sqrt(E_h[0]**2 + E_h[1]**2) * 1e-4
                dQ_h *= np.exp(-dt * step / t_h)

            if t_r:
                dQ_h *= np.exp(-dt * step / (t_r / 2.2))

            dQ_h_step[sel_h] += dQ_h

            Q_ind_tot_h[sel_h] += dQ_h

            I_ind_tot[step] += dQ_h.sum() / dt

        # Store
        store_if_needed(step, next_step, Q_ind_tot_e=Q_ind_tot_e,
                        Q_ind_tot_h=Q_ind_tot_h, dt=dt,
                        dQ_e=dQ_e, dQ_h=dQ_h, T=T, I_ind_e=I_ind_e,
                        I_ind_h=I_ind_h, p_e=p_e, p_h=p_h, q_max=q0,
                        i_step=i_step, n_store=n_store,
                        sel_e=sel_e, sel_h=sel_h)

        # Position change in um
        d_p_e, d_p_h = v_e * dt * 1e-5, v_h * dt * 1e-5

        # Update position
        p_e[:, sel_e] = p_e[:, sel_e] + d_p_e
        p_h[:, sel_h] = p_h[:, sel_h] + d_p_h

        # Correct boundaries (e.g. leaving sensor due to diffusion)
        _correct_boundary(geom_descr, pot_descr,
                          x=p_e[0, :], y=p_e[1, :],
                          sel=sel_e, is_electron=True)
        _correct_boundary(geom_descr, pot_descr,
                          x=p_h[0, :], y=p_h[1, :],
                          sel=sel_h, is_electron=False)

        # Check boundaries and update selection
        sel_e, new_e = _in_boundary(geom_descr, pot_descr=pot_descr,
                                    x=p_e[0, :], y=p_e[1, :],
                                    sel=sel_e)
        sel_h, new_h = _in_boundary(geom_descr, pot_descr=pot_descr,
                                    x=p_h[0, :], y=p_h[1, :],
                                    sel=sel_h)

        # Force a storing step at for e-h pairs where one is finished
        next_step[new_e] = step + 1
        next_step[new_h] = step + 1

        p_e[:, ~sel_e] = np.nan
        p_h[:, ~sel_h] = np.nan

        progress_bar.update(step)
    progress_bar.finish()

    return traj_e, traj_h, I_ind_e, I_ind_h, T, I_ind_tot, Q_ind_tot_e, Q_ind_tot_h
