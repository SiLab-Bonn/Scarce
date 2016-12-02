r"""The solver module provides numerical ODE solver functions.

These functions are either implemented here or provide an interface to other solvers.
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
    ''' Solve the drift-diffusion equation for pseudo particles via euler forward difference + Monte Carlo diffusion
    '''
    
    def __init__(self, pot_descr, pot_w_descr,
                 T=300):
        '''
        Parameters
        ----------
        p0 : array_like, shape = (n, 2)
            n positions of e-h-pairs in x, y at t = 0
    
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
    
    def solve(self, p0, dt, n_steps=4000):
        ''' dt in ns
        '''

        # Start positions
        p_e, p_h = p0, p0
        
        dx = 0.1
        E_tot = 0
        
        for x in np.arange(100., 0., -dx):
            E_tot += self.pot_descr.get_field(0., x)[1] * dx
            
        print 'E_tot', E_tot
#         print np.sum([self.pot_descr.get_field(0., x) * dx for x ])

        print p_e.shape
        traj = np.zeros(shape=(n_steps, p_e.shape[0], p_e.shape[1]))
        
        progress_bar = progressbar.ProgressBar(widgets=['', progressbar.Percentage(), ' ', 
                                                        progressbar.Bar(marker='*', left='|', right='|'), ' ', 
                                                        progressbar.AdaptiveETA()], 
                                               maxval=n_steps, term_width=80)
        progress_bar.start()
        
        for step in range(n_steps):                
            # Electric field in V/cm
            E_e = self.pot_descr.get_field(p_e[0], p_e[1]) * 1e4
            E_h = self.pot_descr.get_field(p_h[0], p_h[1]) * 1e4

            # Mobility in cm2 / Vs
            mu_e = silicon.get_mobility(np.sqrt(E_e[0] ** 2 + E_e[1] ** 2), temperature=self.T, is_electron=True)
            mu_h = silicon.get_mobility(np.sqrt(E_h[0] ** 2 + E_h[1] ** 2), temperature=self.T, is_electron=False)

            # Velocity in cm / s
            v_e, v_h = -E_e * mu_e, E_h * mu_h
            
            # Position change in um
            d_p_e, d_p_h = v_e * dt * 1e-5, v_h * dt * 1e-5

            # Update position
            p_e, p_h = p_e + d_p_e, p_h + d_p_h
            
            traj[step] = p_e
            
            progress_bar.update(step)
        progress_bar.finish()
            
        print traj.shape 
        plt.plot(np.arange(n_steps) * dt, traj[:, 1, 0])
        plt.xlabel('t [ns]')
        plt.show()