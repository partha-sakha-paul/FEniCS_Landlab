from dolfinx.nls.petsc import NewtonSolver
from dolfinx.fem.petsc import NonlinearProblem

def setup_solver(weak_form, W, W_prev, V_test, bcs):
    """
    Setup a Newton solver for solving the nonlinear shallow water equations.

    Parameters:
    - weak_form (callable): Function returning the weak form of the problem.
    - W (Function): Current solution function.
    - W_prev (Function): Previous time-step solution function.
    - V_test (FunctionSpace): Test function space.
    - bcs (list): List of boundary conditions.

    Returns:
    - solver (NewtonSolver): Configured Newton solver.
    """

    # Define the nonlinear problem using the weak form
    problem = NonlinearProblem(weak_form(W, W_prev, V_test, wd=True), W, bcs=bcs)

    # Create and configure the Newton solver
    solver = NewtonSolver(W.function_space.mesh.comm, problem)
    solver.rtol = 1e-6  # Relative tolerance
    solver.atol = 1e-6  # Absolute tolerance
    solver.max_it = 50  # Maximum number of iterations
    solver.error_on_nonconvergence = True  # Raise an error if solver does not converge

    return solver