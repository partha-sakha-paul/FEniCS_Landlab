import numpy as np
from dolfinx import mesh, fem
from petsc4py.PETSc import ScalarType

def open_boundary(x):
    return np.isclose(x[0], 0)

def wall_boundary(x, Ly):
    return np.logical_not(np.isclose(x[0], 0)) | np.isclose(x[1], 0) | np.isclose(x[1], Ly)

def setup_boundary_conditions(domain, V, Ly):
    boundary_facets_open = mesh.locate_entities_boundary(domain, domain.topology.dim - 1, open_boundary)
    boundary_facets_wall = mesh.locate_entities_boundary(domain, domain.topology.dim - 1, lambda x: wall_boundary(x, Ly))

    tidal_value = fem.Constant(domain, ScalarType(0.0))
    sub_dofs_h_left = fem.locate_dofs_topological(V.sub(0), domain.topology.dim - 1, boundary_facets_open)
    bc_h_left = fem.dirichletbc(tidal_value, sub_dofs_h_left, V.sub(0))

    sub_dofs_ux_others = fem.locate_dofs_topological(V.sub(1), domain.topology.dim - 1, boundary_facets_wall)
    sub_dofs_uy_others = fem.locate_dofs_topological(V.sub(2), domain.topology.dim - 1, boundary_facets_wall)
    bc_ux_others = fem.dirichletbc(fem.Constant(domain, ScalarType(0.0)), sub_dofs_ux_others, V.sub(1))
    bc_uy_others = fem.dirichletbc(fem.Constant(domain, ScalarType(0.0)), sub_dofs_uy_others, V.sub(2))

    return [bc_h_left, bc_ux_others, bc_uy_others], tidal_value
