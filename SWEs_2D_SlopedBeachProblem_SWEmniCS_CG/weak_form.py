import ufl
from ufl import TrialFunction, TestFunction, dx, inner
from ufl import avg, jump, dot, dS, as_vector, conditional, sqrt, as_tensor, grad
from dolfinx import fem
from petsc4py.PETSc import ScalarType
import numpy as np
from dolfinx.fem.petsc import assemble_vector
from dolfinx.fem import form

def wd_correction(h, alpha_sq):
    """
    Compute correction to water column height for wetting-and-drying stabilization.
    
    Parameters:
    h : ufl.Expr
        Water depth variable.
    alpha_sq : float
        Squared wetting-drying parameter.

    Returns:
    ufl.Expr
        Corrected water depth.
    """
    return 0.5 * (ufl.sqrt(h**2 + alpha_sq) - h)
    
def get_standard_vars(u, z, wd_alpha_sq, form='H', wd=True):
    """Return a standardized representation of the solution variables.

    Parameters:
    u (tuple): A tuple containing (H, ux, uy), where:
        - H (float): Water height or total head.
        - ux (float): Velocity in the x-direction.
        - uy (float): Velocity in the y-direction.
    z (float): Bottom elevation.
    wd_alpha_sq (float): Correction factor for wetting and drying.
    form (str, optional): The desired output form. Choices:
        - 'H' (default): Returns (H, ux, uy).
        - 'h': Returns (h, ux, uy), where h is water depth.
        - 'flux': Returns (H, Hux, Huy).
    wd (bool, optional): Whether to apply the wetting and drying correction.

    Returns:
    tuple: The standardized variables based on the chosen form.

    Raises:
    ValueError: If an invalid output form is provided.
    """
    H, ux, uy = u[0], u[1], u[2]  # Unpacking input variables
    h = H - z  # Compute water depth

    if wd:
        if form == 'H' or form == 'flux':  # Apply correction for wetting and drying if applicable
            H = H + wd_correction(H, wd_alpha_sq)
    else:
        print("WD NONACTIVE")  # Debugging statement for when wetting/drying is disabled

    Hux, Huy = H * ux, H * uy  # Compute momentum fluxes

    # Return the requested form of variables
    if form == 'H': 
        return H, ux, uy
    elif form == 'h': 
        return h, ux, uy
    elif form == 'flux': 
        return H, Hux, Huy
    else:
        raise ValueError(f"Invalid output form '{form}'")  # Handle incorrect input
    
def make_Fu(u, z, wd_alpha_sq, g, wd=True):
    """Compute the flux tensor for the shallow water equations.

    Parameters:
    u (tuple): A tuple containing (H, ux, uy), where:
        - H (float): Water height or total head.
        - ux (float): Velocity in the x-direction.
        - uy (float): Velocity in the y-direction.
    z (float): Bottom elevation.
    wd_alpha_sq (float): Correction factor for wetting and drying.
    g (float): Gravitational acceleration.
    wd (bool, optional): Whether to apply the wetting and drying correction (default: True).

    Returns:
    Tensor: A 3x2 tensor representing the flux terms.
    """
    # Get standardized representations of variables
    H, ux, uy = get_standard_vars(u, z, wd_alpha_sq, form='H', wd=True)
    h, _, _ = get_standard_vars(u, z, wd_alpha_sq, form='h', wd=True)

    # Compute corrected bottom elevation
    z_c = z + wd_correction(u[0], wd_alpha_sq)

    # Construct the flux tensor using the computed values
    return as_tensor([
        [(H) * ux, (H) * uy],  # Mass flux
        [(H) * ux * ux + 0.5 * g * (h + 2 * z_c) * h, (H) * ux * uy],  # Momentum flux in x-direction
        [(H) * ux * uy, (H) * uy * uy + 0.5 * g * (h + 2 * z_c) * h]   # Momentum flux in y-direction
    ])

def make_Fu_wall(u, z, wd_alpha_sq, g, wd=True):
    """Compute the flux tensor for the shallow water equations.

    Parameters:
    u (tuple): A tuple containing (H, ux, uy), where:
        - H (float): Water height or total head.
        - ux (float): Velocity in the x-direction.
        - uy (float): Velocity in the y-direction.
    z (float): Bottom elevation.
    wd_alpha_sq (float): Correction factor for wetting and drying.
    g (float): Gravitational acceleration.
    wd (bool, optional): Whether to apply the wetting and drying correction (default: True).

    Returns:
    Tensor: A 3x2 tensor representing the flux terms.
    """
    # Get standardized representations of variables
    H, ux, uy = get_standard_vars(u, z, wd_alpha_sq, form='H', wd=True)
    h, _, _ = get_standard_vars(u, z, wd_alpha_sq, form='h', wd=True)

    # Compute corrected bottom elevation
    z_c = z + wd_correction(u[0], wd_alpha_sq)

    # Construct the flux tensor using the computed values
    return as_tensor([
        [0, 0],  # Mass flux
        [0.5 * g * (h + 2 * z_c) * h, 0],  # Momentum flux in x-direction
        [0, 0.5 * g * (h + 2 * z_c) * h]   # Momentum flux in y-direction
    ])

def get_friction(u, z, g, wd_alpha_sq, wd=True):
    """Compute the friction term for the shallow water equations.

    Parameters:
    u (tuple): A tuple containing (H, ux, uy), where:
        - H (float): Water height or total head.
        - ux (float): Velocity in the x-direction.
        - uy (float): Velocity in the y-direction.
    z (float): Bottom elevation.
    g (float): Gravitational acceleration.
    wd_alpha_sq (float): Correction factor for wetting and drying.
    wd (bool, optional): Whether to apply the wetting and drying correction (default: True).

    Returns:
    as_vector: A 3-element vector representing frictional forces in the system.
    """
    # Get standardized variables, ensuring H is correctly adjusted
    H, ux, uy = get_standard_vars(u, z, wd_alpha_sq, form='H', wd=True)

    # Small threshold to avoid division by zero or undefined behavior
    eps = 1e-8
    mannings = 0.02  # Manning's roughness coefficient

    # Compute velocity magnitude, ensuring it is nonzero when below the threshold
    mag_v = conditional(pow(ux * ux + uy * uy, 0.5) < eps, 0, pow(ux * ux + uy * uy, 0.5))

    # Compute and return frictional forces
    return as_vector((
        0,  # No frictional force in mass balance
        g * mannings * mannings * ux * mag_v * pow(H, -1 / 3),  # Friction in x-direction
        g * mannings * mannings * uy * mag_v * pow(H, -1 / 3)   # Friction in y-direction
    ))
    
def make_Source(u, z, wd_alpha_sq, g, wd=True):
    """Compute the source term for the shallow water equations.

    Parameters:
    u (tuple): A tuple containing (H, ux, uy), where:
        - H (float): Water height or total head.
        - ux (float): Velocity in the x-direction.
        - uy (float): Velocity in the y-direction.
    z (float): Bottom elevation.
    wd_alpha_sq (float): Correction factor for wetting and drying.
    g (float): Gravitational acceleration.
    wd (bool, optional): Whether to apply the wetting and drying correction (default: True).

    Returns:
    as_vector: A 3-element vector representing the source terms due to gravity and friction.
    """
    # Get standardized variables
    H, ux, uy = get_standard_vars(u, z, wd_alpha_sq, form='H', wd=True)

    # Compute corrected bottom elevation
    z_c = z + wd_correction(u[0], wd_alpha_sq)

    # Compute gravitational source term if wetting and drying is active
    if wd:
        h, _, _ = get_standard_vars(u, z, wd_alpha_sq, form='h', wd=True)
        g_vec = as_vector((
            0,                # No source term for mass conservation
            -g * h * z_c.dx(0),  # Gravity effect in x-direction
            -g * h * z_c.dx(1)   # Gravity effect in y-direction
        ))
    else:
        g_vec = as_vector((0, 0, 0))  # No gravity effect if wetting/drying is inactive

    # Compute total source term by adding friction effects
    source = g_vec + get_friction(u, z, g, wd_alpha_sq, wd=True)

    return source

def add_bcs_to_weak_form(u, W_open, z, g, wd_alpha_sq, domain, ds, F_u, F_u_open, F_u_wall, v, F, wd=True):
    """Apply boundary conditions to the weak form of the shallow water equations.

    Parameters:
    u (tuple): Solution variables (H, ux, uy).
    W_open (tuple): External state variables at open boundary.
    z (float): Bottom elevation.
    g (float): Gravitational acceleration.
    wd_alpha_sq (float): Wetting and drying correction factor.
    domain: Computational domain for the problem.
    ds: Boundary integration measure.
    F_u: Flux tensor computed using the current state variables.
    F_u_open: Flux tensor computed using the external state at open boundary.
    F_u_wall: Flux tensor computed using the external state at closed boundary.
    v: Test function.
    F: Weak form equation that will be modified.
    wd (bool, optional): Whether to apply wetting and drying correction (default: True).

    Returns:
    Updated weak form equation F with applied boundary conditions.
    """
    # Compute unit normal to the domain boundary
    n = ufl.FacetNormal(domain)

    # Apply open boundary conditions to the weak form
    F += dot(dot(F_u_open, n), v) * ds(2)

    # Apply wall boundary conditions to the weak form
    F += dot(dot(F_u_wall, n), v) * ds(1)

    return F

def weak_form(U, U_prev, U_prevprev, V_test, W_open, domain, z, dt, ds, g, mannings, step, wd_alpha=0.36, theta1=1.0, wd=True):
    """
    Constructs the weak form of the shallow water equations for time-stepping.

    Parameters:
    U (tuple): Current state variables (H, ux, uy).
    U_prev (tuple): State variables from the previous time step.
    U_prevprev (tuple): State variables from two time steps ago.
    V_test (tuple): Test functions (v_h, v_ux, v_uy).
    W_open (tuple): Open boundary state variables (H_open, ux_open, uy_open).
    domain: Computational domain.
    z (float): Bottom elevation.
    dt (float): Time step size.
    ds: Surface measure for boundary integrals.
    g (float): Gravitational acceleration.
    mannings (float): Manning's roughness coefficient.
    step (int): Current time step number.
    wd_alpha (float, optional): Wetting and drying correction factor (default: 0.36).
    theta1 (float, optional): Theta scheme parameter (default: 1.0).
    wd (bool, optional): Whether to apply wetting and drying correction (default: True).

    Returns:
    F: The weak form equation to be solved.
    """

    # Extract the primary variables from the solution states
    H, ux, uy = U
    H_prev, ux_prev, uy_prev = U_prev
    H_prevprev, ux_prevprev, uy_prevprev = U_prevprev
    H_open, ux_open, uy_open = W_open
    v_h, v_ux, v_uy = V_test  # Test functions for weak formulation

    # Compute squared wetting and drying correction factor
    wd_alpha_sq = ScalarType(wd_alpha**2)

    # Adjust theta1 for the first two steps (explicit for step < 2, implicit otherwise)
    # theta1 = 1.0 if step < 2 else 1.0
    theta1 = fem.Constant(domain, ScalarType(theta1))

    # Compute unit normal vector at element boundaries
    n = ufl.FacetNormal(domain)

    # Define test function as a vector
    v = as_vector((v_h, v_ux, v_uy))

    # Compute flux variables for current and previous states
    Q = as_vector(get_standard_vars(U, z, wd_alpha_sq, form='flux', wd=True))
    Q_prev = as_vector(get_standard_vars(U_prev, z, wd_alpha_sq, form='flux', wd=True))
    Q_prevprev = as_vector(get_standard_vars(U_prevprev, z, wd_alpha_sq, form='flux', wd=True))

    # Compute time derivative using a second-order backward difference formula (BDF2)
    dQdt = (
        theta1 * fem.Constant(domain, ScalarType(1/dt)) * (1.5 * Q - 2 * Q_prev + 0.5 * Q_prevprev) + 
        (1 - theta1) * fem.Constant(domain, ScalarType(1/dt)) * (Q - Q_prev)
    )

    # Debugging print statement to check solution array
    # print('U for Fu', U.sub(0).x.array)
    u_open = as_vector((H_open,U[1],U[2]))

    # Compute the flux tensor for the state variables
    F_u = make_Fu(U, z, wd_alpha_sq, g, wd=True)
    F_u_open = make_Fu(u_open, z, wd_alpha_sq, g, wd=True)
    F_u_wall = make_Fu_wall(U, z, wd_alpha_sq, g, wd=True)
    # a = -inner(F_u_open,grad(v))*dx
    # print('fuopen norm', np.linalg.norm(assemble_vector(form(a)).array))
    # Compute the source term (gravitational and frictional forces)
    S = make_Source(U, z, wd_alpha_sq, g, wd=True)

    # Initialize weak form with the divergence of the flux tensor
    F = -inner(F_u, grad(v)) * dx
    # print('weakF1', np.linalg.norm(assemble_vector(form(F)).array))

    # Apply boundary conditions to the weak form
    F = add_bcs_to_weak_form(U, W_open, z, g, wd_alpha_sq, domain, ds, F_u, F_u_open, F_u_wall, v, F, wd=True)
    # print('weakF2', np.linalg.norm(assemble_vector(form(F)).array))
    
    # Add source term contribution
    F += inner(S, v) * dx
    # print('weakF3', np.linalg.norm(assemble_vector(form(F)).array))
    
    # Add time derivative term
    F += inner(dQdt, v) * dx
    # print('weakF4', np.linalg.norm(assemble_vector(form(F)).array))
    
    # Print norm of the assembled weak form for debugging at the first time step
    if step == 0:
        print('weakF norm', np.linalg.norm(assemble_vector(form(F)).array))

    return F