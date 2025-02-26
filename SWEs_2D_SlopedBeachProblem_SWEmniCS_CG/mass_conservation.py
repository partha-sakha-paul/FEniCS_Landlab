import numpy as np

def mass_conservation_residual(H, u, H_prev, u_prev, dx, dy, dt, maps):
    """
    Computes the mass conservation residual for the shallow water equations with a mixed DG(1) element.

    Parameters:
    - H: Current time step water depth (dolfinx.fem.Function).
    - u: Current time step velocity vector (dolfinx.fem.Function, shape=(2,)).
    - H_prev: Previous time step water depth.
    - u_prev: Previous time step velocity vector.
    - dx, dy: Spatial step sizes in x and y directions.
    - dt: Time step size.
    - maps: Index mappings for extracting DOF values from functions.

    Returns:
    - residual: Numpy array of residual values representing the mass conservation error.
    """

    # Extract values of H and H_prev using maps[0]
    H_values = H.x.array[maps[0]]
    H_prev_values = H_prev.x.array[maps[0]]

    # Extract velocity components from the vector field ux and uy separately using maps[1]
    ux_values = u.x.array[maps[1][::2]]  # Every first element (x-component)
    uy_values = u.x.array[maps[1][1::2]]  # Every second element (y-component)

    # Initialize residual array
    residual = np.zeros_like(H_values)

    # Compute time derivative term
    residual += (H_values - H_prev_values) / dt

    # Compute spatial derivative terms using finite difference approximations
    residual += np.gradient(H_values * ux_values, dx)
    residual += np.gradient(H_values * uy_values, dy)

    return residual