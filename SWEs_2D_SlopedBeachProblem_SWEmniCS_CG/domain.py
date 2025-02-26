import numpy as np
from mpi4py import MPI
from dolfinx import mesh, fem
import basix

def create_domain(Lx, Ly, nx, ny):
    """
    Create a 2D triangular mesh for a rectangular domain.

    Parameters:
    - Lx, Ly: Length of the domain in x and y directions.
    - nx, ny: Number of elements in x and y directions.

    Returns:
    - domain: The created DOLFINx mesh object.
    """
    domain = mesh.create_rectangle(
        MPI.COMM_WORLD,  # Use MPI for parallel computing
        [np.array([0.0, 0.0]), np.array([Lx, Ly])],  # Bottom-left and top-right corners
        [nx, ny],  # Number of divisions in x and y directions
        cell_type=mesh.CellType.triangle  # Use triangular elements
    )
    
    # Create topological connectivity between mesh entities
    domain.topology.create_connectivity(domain.topology.dim - 1, domain.topology.dim)  # Facets to cells
    domain.topology.create_connectivity(0, 1)  # Vertices to edges

    return domain

def create_function_space(domain):
    """
    Create a mixed finite element function space on the given mesh.

    Parameters:
    - domain: The DOLFINx mesh object.

    Returns:
    - Function space consisting of three continuous Galerkin (CG) elements of degree 1.
    """
    element = basix.ufl.mixed_element(
        [basix.ufl.element("CG", str(domain.ufl_cell()), 1), basix.ufl.element("CG", str(domain.ufl_cell()), 1, shape=(2,))]
    )

    # Create the function space on the given mesh
    V = fem.functionspace(domain, element)
    return V
