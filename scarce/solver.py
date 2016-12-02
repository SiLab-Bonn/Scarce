r"""The solver module provides numerical ODE solver functions.

These functions are either implemented here or provide an interface to other solvers.
"""

import fipy

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


def solve_drift_diffusion(x, y, pot_descr, w_pot_descr):
    ''' Solve the drift-diffusion equation for pseudo particles via euler forward difference + Monte Carlo diffusion.

    Parameters
    ----------
    x : array_like
        x positions of the pseudo e, h pairs

    y : array_like
        y positions of the pseudo e, h pairs

    pot_descr : fields.Description
        Description object to describe the potential.

    w_pot_descr : fields.Description
        Description object to describe the weighting potential.

    Notes
    -----
    
    '''
    pass