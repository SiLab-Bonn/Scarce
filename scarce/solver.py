r"""The solver module takes care of the correct settings for the backend solver package used.

This is important to ensure correct solutions on different operating systems.

"""

import fipy


def solve(var, equation, **kwargs):   
    # Force a linear solver
    # The others do not always converge (?!)
    solver = fipy.solvers.LinearLUSolver

    # Reasonable and more strict convergence criteria
    tolerance = 1e-15
    iterations = 10000

    equation.solve(var=var,
                   solver=solver(tolerance=tolerance, iterations=iterations),
                   **kwargs)



