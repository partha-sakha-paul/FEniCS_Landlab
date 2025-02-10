import numpy as np

def mass_conservation_residual(h, ux, uy, h_prev, ux_prev, uy_prev, dx, dy, dt, maps):
    """
    Computes the mass conservation residual for the shallow water equations.

    Parameters:
    - h, ux, uy: Current time step water depth and velocity components (dolfinx.fem.Function).
    - h_prev, ux_prev, uy_prev: Previous time step water depth and velocity components.
    - dx, dy: Spatial step sizes in x and y directions.
    - dt: Time step size.
    - maps: Index mappings for extracting DOF values from functions.

    Returns:
    - residual: Numpy array of residual values representing the mass conservation error.

    This function computes the residual of the mass conservation equation
    using finite difference approximations for the derivatives.
    """

    # Extract values of h, ux, uy at the current and previous time steps
    h_values = h.x.array[maps[0]]
    ux_values = ux.x.array[maps[1]]
    uy_values = uy.x.array[maps[2]]

    h_prev_values = h_prev.x.array[maps[0]]
    ux_prev_values = ux_prev.x.array[maps[1]]
    uy_prev_values = uy_prev.x.array[maps[2]]

    # Initialize residual array
    residual = np.zeros_like(h_values)

    # Compute time derivative term
    residual += (h_values - h_prev_values) / dt

    # Compute spatial derivative terms using finite difference approximations
    residual += np.gradient(h_values * ux_values, dx)  
    residual += np.gradient(h_values * uy_values, dy) 

    return residual
