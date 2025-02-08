import numpy as np

def mass_conservation_residual(h, ux, uy, h_prev, ux_prev, uy_prev, dx, dy, dt, maps):
    """
    Compute the residual of the mass conservation equation.

    Parameters:
    - h: Current water depth field.
    - ux: Current x-component of velocity.
    - uy: Current y-component of velocity.
    - h_prev: Previous time step water depth.
    - ux_prev: Previous time step x-component of velocity.
    - uy_prev: Previous time step y-component of velocity.
    - dx: Grid spacing in the x-direction.
    - dy: Grid spacing in the y-direction.
    - dt: Time step size.
    - maps: List of index mappings for extracting field values.

    Returns:
    - residual: Residual of the mass conservation equation.
    """

    # Extract values from function arrays based on given maps
    h_values = h.x.array[maps[0]]
    ux_values = ux.x.array[maps[1]]
    uy_values = uy.x.array[maps[2]]

    h_prev_values = h_prev.x.array[maps[0]]
    ux_prev_values = ux_prev.x.array[maps[1]]
    uy_prev_values = uy_prev.x.array[maps[2]]

    # Initialize residual array with the same shape as h_values
    residual = np.zeros_like(h_values)

    # Compute temporal change in water depth
    residual += (h_values - h_prev_values) / dt

    # Compute spatial derivatives using numerical gradients
    residual += np.gradient(h_values * ux_values, dx)  
    residual += np.gradient(h_values * uy_values, dy)  

    return residual
