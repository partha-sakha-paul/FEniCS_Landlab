import numpy as np
from dolfinx import mesh, fem
from petsc4py.PETSc import ScalarType

def open_boundary(x):
    """
    Function to identify the open boundary (left edge at x = 0).

    Parameters:
    - x: Array of coordinates.

    Returns:
    - Boolean array indicating points that lie on the open boundary.
    """
    return np.isclose(x[0], 0)

def wall_boundary(x, Ly):
    """
    Function to identify the wall boundaries, which include:
    - The top and bottom boundaries (y = 0 and y = Ly).
    - The right boundary (x ≠ 0).
    
    Parameters:
    - x: Array of coordinates.
    - Ly: Height of the domain.

    Returns:
    - Boolean array indicating points that lie on the wall boundaries.
    """
    return np.logical_not(np.isclose(x[0], 0)) | np.isclose(x[1], 0) | np.isclose(x[1], Ly)

def setup_boundary_conditions(domain, V, Ly):
    """
    Set up Dirichlet boundary conditions for the problem.

    Parameters:
    - domain: The computational mesh.
    - V: The mixed function space (containing h, ux, uy).
    - Ly: Height of the domain.

    Returns:
    - List of Dirichlet boundary conditions applied to h, ux, and uy.
    - The tidal value as a fem.Constant.
    """
    # Locate boundary facets for open boundary (x = 0)
    boundary_facets_open = mesh.locate_entities_boundary(domain, domain.topology.dim - 1, open_boundary)

    # Locate boundary facets for wall boundaries (y = 0, y = Ly, or x ≠ 0)
    boundary_facets_wall = mesh.locate_entities_boundary(domain, domain.topology.dim - 1, lambda x: wall_boundary(x, Ly))

    # Define tidal boundary value (set to 0 for now, can be updated dynamically)
    tidal_value = fem.Constant(domain, ScalarType(0.0))

    # Apply Dirichlet BC for h at the open boundary (x = 0)
    sub_dofs_h_left = fem.locate_dofs_topological(V.sub(0), domain.topology.dim - 1, boundary_facets_open)
    bc_h_left = fem.dirichletbc(tidal_value, sub_dofs_h_left, V.sub(0))

    # Apply Dirichlet BC for ux and uy at the wall boundaries (no-slip condition)
    sub_dofs_ux_others = fem.locate_dofs_topological(V.sub(1), domain.topology.dim - 1, boundary_facets_wall)
    sub_dofs_uy_others = fem.locate_dofs_topological(V.sub(2), domain.topology.dim - 1, boundary_facets_wall)

    bc_ux_others = fem.dirichletbc(fem.Constant(domain, ScalarType(0.0)), sub_dofs_ux_others, V.sub(1))
    bc_uy_others = fem.dirichletbc(fem.Constant(domain, ScalarType(0.0)), sub_dofs_uy_others, V.sub(2))

    # Return the list of boundary conditions and the tidal value for further modifications
    return [bc_h_left, bc_ux_others, bc_uy_others], tidal_value
