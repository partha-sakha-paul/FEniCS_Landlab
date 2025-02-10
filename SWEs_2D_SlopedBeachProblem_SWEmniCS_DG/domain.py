import numpy as np
from mpi4py import MPI
from dolfinx import mesh, fem
import basix

def create_domain(Lx=13800.0, Ly=7200.0, nx=12, ny=6):
    """
    Create a 2D rectangular mesh using triangular elements.

    Parameters:
    Lx, Ly : float
        Length of the domain in the x and y directions.
    nx, ny : int
        Number of divisions (cells) along x and y axes.

    Returns:
    domain : dolfinx.mesh.Mesh
        The generated mesh.
    """
    domain = mesh.create_rectangle(
        MPI.COMM_WORLD,  # Use MPI for parallel processing
        [np.array([0.0, 0.0]), np.array([Lx, Ly])],  # Define lower-left and upper-right corners
        [nx, ny],  # Number of cells in x and y directions
        cell_type=mesh.CellType.triangle  # Use triangular elements
    )

    # Create connectivity between mesh entities (faces, edges, and vertices)
    domain.topology.create_connectivity(domain.topology.dim - 1, domain.topology.dim)  # Face-to-cell
    domain.topology.create_connectivity(0, 1)  # Vertex-to-edge

    return domain

def create_function_space(domain):
    """
    Create a mixed finite element function space with three components.

    Parameters:
    domain : dolfinx.mesh.Mesh
        The computational domain.

    Returns:
    V : dolfinx.fem.FunctionSpace
        The mixed function space.
    """
    # Define a mixed finite element with three DG(1) components
    element = basix.ufl.mixed_element(
        [basix.ufl.element("DG", str(domain.ufl_cell()), 1)] * 3
    )

    # Create the function space on the given mesh
    V = fem.functionspace(domain, element)

    return V
