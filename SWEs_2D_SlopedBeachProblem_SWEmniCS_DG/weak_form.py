import ufl
from ufl import TrialFunction, TestFunction, dx, inner
from ufl import avg, jump, dot, dS, as_vector, conditional, sqrt, as_tensor
from dolfinx import fem
from petsc4py.PETSc import ScalarType

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

def weak_form(U, U_prev, U_prevprev, V_test, W_open, domain, z, dt, ds, g, mannings, step, wd_alpha=0.36, theta1=1.0, wd=True):
    """
    Constructs the weak form of the shallow water equations using a variational formulation.

    Parameters:
    U : tuple (h, ux, uy)
        Current solution variables: water depth and velocity components.
    U_prev : tuple (h_prev, ux_prev, uy_prev)
        Solution variables at the previous time step.
    U_prevprev : tuple (h_prevprev, ux_prevprev, uy_prevprev)
        Solution variables at the second-last time step.
    V_test : tuple (v_h, v_ux, v_uy)
        Test functions for the variational formulation.
    W_open : tuple (h_open, ux_open, uy_open)
        Open boundary condition values.
    domain : dolfinx.mesh.Mesh
        Computational domain.
    z : ufl.Expr
        Bathymetry function (bed elevation).
    dt : float
        Time step size.
    ds : ufl.Measure
        Boundary integration measure.
    g : float
        Acceleration due to gravity.
    mannings : float
        Manning's roughness coefficient.
    step : int
        Current time step index.
    wd_alpha : float, optional
        Wetting-drying parameter (default: 0.36).
    theta1 : float, optional
        Time integration weight (default: 1.0).
    wd : bool, optional
        Whether to apply wetting-drying correction (default: True).

    Returns:
    ufl.Form
        The weak form of the shallow water equations.
    """

    # Extract solution variables
    h, ux, uy = U
    h_prev, ux_prev, uy_prev = U_prev
    h_prevprev, ux_prevprev, uy_prevprev = U_prevprev 
    h_open, ux_open, uy_open = W_open
    v_h, v_ux, v_uy = V_test

    # Apply wetting-drying correction if enabled
    if wd:
        wd_alpha_sq = ScalarType(wd_alpha**2)
        h_corrected = h + wd_correction(h, wd_alpha_sq)
        h_open_corrected = h_open + wd_correction(h_open, wd_alpha_sq)
    else:
        h_corrected = h
        h_open_corrected = h_open

    # Adjust theta1 for startup time steps
    theta1 = 0.0 if step < 2 else 1.0

    # Compute the unit normal vector on the boundary
    n = ufl.FacetNormal(domain)  

    # Time derivatives using backward difference formulas
    dHdt = theta1 * (1.5 * h - 2.0 * h_prev + 0.5 * h_prevprev) / dt + (1 - theta1) * (h - h_prev) / dt
    dHUxdt = theta1 * (1.5 * (h + z) * ux - 2.0 * (h_prev + z) * ux_prev + 0.5 * (h_prevprev + z) * ux_prevprev) / dt + \
             (1 - theta1) * ((h + z) * ux - (h_prev + z) * ux_prev) / dt
    dHUydt = theta1 * (1.5 * (h + z) * uy - 2.0 * (h_prev + z) * uy_prev + 0.5 * (h_prevprev + z) * uy_prevprev) / dt + \
             (1 - theta1) * ((h + z) * uy - (h_prev + z) * uy_prev) / dt

    # Weak form of the continuity equation (conservation of mass)
    continuity = dHdt * v_h * dx - \
                ((h + z) * ux * v_h.dx(0) + (h + z) * uy * v_h.dx(1)) * dx

    # Weak form of the momentum equation in the x-direction
    momentum_x = dHUxdt * v_ux * dx - \
            (((h + z) * ux**2 + 0.5 * g * (h + 2 * z) * h) * v_ux.dx(0) + ((h + z) * ux * uy) * v_ux.dx(1)) * dx - \
             g * h * z.dx(0) * v_ux * dx + \
             g * mannings**2 * (abs(ux) * ux / h_corrected**(1/3)) * v_ux * dx

    # Weak form of the momentum equation in the y-direction
    momentum_y = dHUydt * v_uy * dx - \
            (((h + z) * ux * uy) * v_uy.dx(0) + ((h + z) * uy**2 + 0.5 * g * (h + 2 * z) * h) * v_uy.dx(1)) * dx - \
            g * h * z.dx(1) * v_uy * dx + \
            g * mannings**2 * (abs(uy) * uy / h_corrected**(1/3)) * v_uy * dx

    # The flux tensors F_u represent the physical fluxes in the equations for interior of the domain
    F_u = as_tensor([
    [(h+z) * ux, (h+z) * uy], 
    [(h+z) * ux * ux + 0.5 * g * (h+2*z)*h, (h+z) * ux * uy],
    [(h+z) * ux * uy, (h+z) * uy * uy + 0.5 * g * (h+2*z)*h]
    ])

    # The flux tensors F_u_open represent the physical fluxes in the equations for open boundaries
    F_u_open = as_tensor([
    [(h_open+z) * ux_open, (h_open+z) * uy_open], 
    [(h_open+z) * ux_open * ux_open + 0.5 * g * (h_open+2*z)*h_open, (h_open+z) * ux_open * uy_open],
    [(h_open+z) * ux_open * uy_open, (h_open+z) * uy_open * uy_open + 0.5 * g * (h_open+2*z)*h_open]
    ])
    
    # add bcs to weak form starts::::
    vel = as_vector((ux, uy))  # velocity vector
    un = dot(vel, n)   # normal component of velocity
    eps = 1e-16  # Small positive value to avoid sqrt of zero in velocity magnitude calculations
    vnorm = conditional(dot(vel, vel) > eps, sqrt(dot(vel, vel)), 0.0)  # norm of velocity vector

    # jump_Q_wall handle wall boundary flux corrections.
    jump_Q_wall = as_vector((0, 2 * (h_corrected+z) * un * n[0], 2 * (h_corrected+z) * un * n[1]))

    # C_wall compute stabilization terms based on flow velocity and depth.
    C_wall = vnorm + sqrt(g * (h_corrected+z))

    # Define the wall boundary condition for the state variables
    u_wall = as_vector((
        U[0],  # Water depth remains the same
        U[1] * n[1] * n[1] - U[1] * n[0] * n[0] - 2 * U[2] * n[0] * n[1],  # Modify x-momentum to enforce no-normal flow
        U[2] * n[0] * n[0] - U[2] * n[1] * n[1] - 2 * U[1] * n[0] * n[1]   # Modify y-momentum similarly
    ))

    # Compute the external flux tensor at the wall using the modified state variables
    Fu_wall_ext = as_tensor([
        [(u_wall[0] + z) * u_wall[1], (u_wall[0] + z) * u_wall[2]],  # First row: h*u terms
        [(u_wall[0] + z) * u_wall[1] * u_wall[1] + 0.5 * g * (u_wall[0] + 2*z) * (u_wall[0] + z), (u_wall[0] + z) * u_wall[1] * u_wall[2]],  # Second row: momentum flux in x-direction
        [(u_wall[0] + z) * u_wall[1] * u_wall[2], (u_wall[0] + z) * u_wall[2] * u_wall[2] + 0.5 * g * (u_wall[0] + 2*z) * (u_wall[0] + z)]  # Third row: momentum flux in y-direction
    ])

    # jump_Q_open handle open boundary flux corrections.
    jump_Q_open = as_vector((h_corrected - h_open_corrected, (h_corrected+z) * ux - (h_open_corrected+z) * ux_open, (h_corrected+z) * uy - (h_open_corrected+z) * uy_open))

    # C_open compute stabilization terms based on flow velocity and depth.
    C_open = vnorm + sqrt(g * conditional(h_open_corrected > h_corrected, h_open_corrected+z, h_corrected+z))

    # boundary_flux applies flux corrections at domain boundaries: 2 is for open, 1 is for wall conditions
    boundary_flux = 0.0
    boundary_flux = dot(0.5 * dot(F_u_open, n) + 0.5 * dot(F_u, n), as_vector([v_h, v_ux, v_uy])) * ds(2) \
                      + dot(0.5 * C_open * jump_Q_open, as_vector([v_h, v_ux, v_uy])) * ds(2)
    boundary_flux += dot(0.5 * dot(F_u, n) + 0.5 * dot(Fu_wall_ext, n), as_vector([v_h, v_ux, v_uy])) * ds(1) \
                      + dot(0.5 * C_wall * jump_Q_wall, as_vector([v_h, v_ux, v_uy])) * ds(1)
    # print(ds(1))
    # add bcs to weak form ends:::::

    # Discontinuous Galerkin (DG) stabilization:::::
    eps = 1e-16  # Small threshold to avoid division by zero in velocity magnitude calculations

    # Define velocity vectors on the '+' and '-' sides of an element interface
    vela = as_vector((ux('+'), uy('+')))  # Velocity from the positive side of the interface
    velb = as_vector((ux('-'), uy('-')))  # Velocity from the negative side of the interface

    # Compute velocity magnitudes, ensuring numerical stability with a small threshold
    vnorma = conditional(dot(vela, vela) > eps, sqrt(dot(vela, vela)), 0.0)  
    vnormb = conditional(dot(velb, velb) > eps, sqrt(dot(velb, velb)), 0.0)

    # Compute the stabilization parameter C using the maximum velocity + gravity wave speed
    # Ensures numerical stability by selecting the largest characteristic speed across the interface
    C = conditional(
        (vnorma + sqrt(g * (h_corrected+z)('+'))) > (vnormb + sqrt(g * (h_corrected+z)('-'))),
        (vnorma + sqrt(g * (h_corrected+z)('+'))),
        (vnormb + sqrt(g * (h_corrected+z)('-')))
    )

    # Compute flux term using the average of F_u across the element interface 
    # and a stabilization term proportional to the jump in solution variables.
    flux = dot(avg(F_u), n('+')) + 0.5 * C * jump(as_vector([
        (h_corrected + z),  # Height component
        (h_corrected + z) * ux,  # Momentum in x-direction
        (h_corrected + z) * uy   # Momentum in y-direction
    ]))

    # Compute the DG flux contribution to the weak form by taking the inner product 
    # of the flux term with the jump in the test functions.
    dg_flux_term = inner(flux, jump(as_vector([v_h, v_ux, v_uy]))) * dS
    # h_func, ux_func, uy_func = U.split()  
    # print(f"h min/max: {h_func.x.array.min()}, {h_func.x.array.max()}")
    # print(f"ux min/max: {ux_func.x.array.min()}, {ux_func.x.array.max()}")
    # print(f"uy min/max: {uy_func.x.array.min()}, {uy_func.x.array.max()}")
    return continuity + momentum_x + momentum_y + boundary_flux + dg_flux_term
