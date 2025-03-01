import numpy as np
from dolfinx import fem, mesh
from dolfinx.mesh import locate_entities_boundary
from ufl import Measure

def initialize_boundary_conditions(V, W_open, domain, Lx, Ly):
    """
    Initializes boundary conditions and measures for weak form integration.

    Parameters:
    V : dolfinx.fem.FunctionSpace
        The function space used in the problem.
    W_open : dolfinx.fem.Function
        The function representing open boundary conditions.
    domain : dolfinx.mesh.Mesh
        The computational mesh.
    Ly : float, optional
        The length of the domain in the y-direction (default: 7200.0).

    Returns:
    ds : ufl.Measure
        Measure for integrating over the domain boundaries.
    sub_dofs_h_open_left : numpy.ndarray
        Degrees of freedom for h at the open boundary (x = 0).
    """

    # Define the open boundary: y = 0
    def open_boundary(x):
        return np.isclose(x[1], 0)

    # Define the wall boundary: all other boundaries except y = 0
    def wall_boundary(x):
        return np.logical_not(np.isclose(x[1], 0)) | np.isclose(x[0], 0) | np.isclose(x[0], Lx)

    # Get the function space associated with W_open
    V_boundary = W_open.function_space

    # Locate degrees of freedom for h at the open boundary (x = 0)
    sub_dofs_h_open_left = fem.locate_dofs_geometrical(
        (V_boundary.sub(0), V_boundary.sub(0).collapse()[0]), open_boundary
    )[0]
    print(sub_dofs_h_open_left)
    # Locate boundary facets based on the defined conditions
    boundary_facets_open = locate_entities_boundary(domain, domain.topology.dim - 1, open_boundary)
    boundary_facets_wall = locate_entities_boundary(domain, domain.topology.dim - 1, wall_boundary)

    # Create facet tags: 
    # - Tag 1: Wall boundary (non-open boundary)
    # - Tag 2: Open boundary (x = 0)
    facet_tags = mesh.meshtags(domain, domain.topology.dim - 1, 
                               np.concatenate([boundary_facets_wall, boundary_facets_open]),
                               np.concatenate([np.full(len(boundary_facets_wall), 1), 
                                               np.full(len(boundary_facets_open), 2)]))

    # Define the measure for boundary integration using facet tags
    ds = Measure("ds", domain=domain, subdomain_data=facet_tags)

    return ds, sub_dofs_h_open_left