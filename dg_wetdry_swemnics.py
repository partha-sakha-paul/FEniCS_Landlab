# for calculating residual
def mass_conservation_residual(h, ux, uy, h_prev, ux_prev, uy_prev, dx, dy, dt, maps):
    h_values = h.x.array[maps[0]]
    ux_values = ux.x.array[maps[1]]
    uy_values = uy.x.array[maps[2]]
    h_prev_values = h_prev.x.array[maps[0]]
    ux_prev_values = ux_prev.x.array[maps[1]]
    uy_prev_values = uy_prev.x.array[maps[2]]
    residual = np.zeros_like(h_values)
    residual += (h_values - h_prev_values) / dt
    residual += np.gradient(h_values * ux_values, dx)  
    residual += np.gradient(h_values * uy_values, dy)  
    return residual
    
def find_values_at_vertices(function, mesh):
    """
    Evaluate a DOLFINx Function at all mesh vertices.

    Parameters:
        function (dolfinx.fem.Function): The function to evaluate (e.g., 'z').
        mesh (dolfinx.mesh.Mesh): The mesh containing the vertices.

    Returns:
        list: A list of values of the function at all vertices.
    """
    # Create connectivity between vertices (entities of dimension 0) and cells (entities of topological dimension)
    mesh.topology.create_connectivity(0, mesh.topology.dim)
    # Create connectivity between vertices (dimension 0) and edges (dimension 1)
    mesh.topology.create_connectivity(0, 1)

    # Get the connectivity information (which cells are connected to each vertex)
    cells = mesh.topology.connectivity(0, mesh.topology.dim)

    # Get the coordinates of all vertices in the mesh
    vertices = mesh.geometry.x

    # Initialize an empty list to store function values at vertices
    values_at_vertices = []

    # Iterate over all vertices in the mesh
    for vertex_index, vertex_coords in enumerate(vertices):
        # Create an array for the current vertex coordinates in 3D (adding z=0.0 for 2D meshes)
        points = np.array([[vertex_coords[0], vertex_coords[1], 0.0]])
        # Find the cells associated with the current vertex
        associated_cells = cells.links(vertex_index)
        # If there are any associated cells (not empty)
        if len(associated_cells) > 0:
            # Evaluate the function at the vertex using the first associated cell
            value = function.eval(points, [associated_cells[0]])
            values_at_vertices.append(value)
    return values_at_vertices

import numpy as np
from mpi4py import MPI
from dolfinx import mesh, fem, nls, io
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from dolfinx.fem import locate_dofs_geometrical
from petsc4py.PETSc import ScalarType
from ufl import TrialFunction, TestFunction, div, dx
from ufl.finiteelement import FiniteElement, MixedElement
import ufl
import basix
from dolfinx.plot import vtk_mesh
import pyvista as pv
import os
import matplotlib.pyplot as plt
from dolfinx.fem import assemble_scalar
from ufl import inner, dx
from dolfinx.fem import form
from dolfinx.mesh import locate_entities_boundary
from ufl import Measure
os.makedirs("frames", exist_ok=True)

# domain and mesh creation
Lx, Ly = 13800.0, 7200.0  
nx, ny = 12, 6           
domain = mesh.create_rectangle(
    MPI.COMM_WORLD, 
    [np.array([0.0, 0.0]), np.array([Lx, Ly])], 
    [nx, ny],
    cell_type=mesh.CellType.triangle
)
domain.topology.create_connectivity(domain.topology.dim - 1, domain.topology.dim)
domain.topology.create_connectivity(0, 1)

# parameters
g = 9.81  
mannings = 0.02
t, dt, T = 0.0, 600, 7*24*3600.0
num_steps = int(T / dt)

# element and functionspace
element = basix.ufl.mixed_element([basix.ufl.element("DG", str(domain.ufl_cell()), 1)] * 3)
V = fem.functionspace(domain, element)

# collapsing the subspaces for later use
num_subs = V.num_sub_spaces
spaces = []
maps = []
for i in range(num_subs):
    space_i, map_i = V.sub(i).collapse()
    spaces.append(space_i)
    maps.append(map_i)

# creating bathymetry
V0, _ = V.sub(0).collapse()  
z = fem.Function(V0) 
z.interpolate(lambda x: 5.0 / 13800 * (13800 - x[0]))
z_values_at_vertices = find_values_at_vertices(z, domain)
print(z_values_at_vertices)
print(len(z_values_at_vertices))
z_values_array = np.concatenate(z_values_at_vertices)
print(z_values_array)
x = domain.geometry.x
z_vertex_values = np.zeros(x.shape[0])  
for i, vertex in enumerate(x):         
    z_vertex_values[i] = 5.0 / 13800 * (13800 - vertex[0])  
x[:, 2] = z_vertex_values
print('xshp',x.shape)
print('z.x.array',z.x.array.shape)
print('x[:,2]',x[:,2])

# current functions
W = fem.Function(V)  
h, ux, uy = W.split()

# initializing current functions
h.interpolate(lambda x: np.full_like(x[0], 0.0))
ux.interpolate(lambda x: np.full_like(x[0], 0.0))
uy.interpolate(lambda x: np.full_like(x[0], 0.0))
print('h shape:', h.x.array.shape)

# previous functions and initialization
W_prev = fem.Function(V)  
h_prev, ux_prev, uy_prev = W_prev.split()
h_prev.interpolate(lambda x: np.full_like(x[0], 0.0))
ux_prev.interpolate(lambda x: np.full_like(x[0], 0.0))
uy_prev.interpolate(lambda x: np.full_like(x[0], 0.0))

# trial and test functions
U = TrialFunction(V)
V_test = TestFunction(V)

# wetting - drying correction for dry regions
def wd_correction(h, alpha_sq):
    """Compute correction to water column height for wetting-and-drying."""
    return 0.5 * (ufl.sqrt(h**2 + alpha_sq) - h)
wd_alpha = 0.36
wd_alpha_sq = ScalarType(wd_alpha**2)

# weak form for dg implemetation
def weak_form(U, U_prev, V_test, wd=True):
    h, ux, uy = U
    h_prev, ux_prev, uy_prev = U_prev
    v_h, v_ux, v_uy = V_test
    if wd:
        h_corrected = h + wd_correction(h, wd_alpha_sq)
    else:
        h_corrected = h
    n = ufl.FacetNormal(domain)  
    continuity = (h - h_prev) / dt * v_h * dx - \
                ((h + z) * ux * v_h.dx(0) + (h + z) * uy * v_h.dx(1)) * dx + \
                ((h + z) * ufl.dot(ufl.as_vector([ux, uy]), n) * v_h) * ds(2) + \
                tidal_value * ufl.dot(n, n) * v_h * ds(2) + \
                ((h + z) * ufl.dot(ufl.as_vector([ux, uy]), n) * v_h) * ds(1)  
    
    momentum_x = ((h+z)*ux - (h_prev+z)*ux_prev) / dt * v_ux * dx - \
             (ux * (h + z) * ux * v_ux.dx(0) + ux * (h + z) * uy * v_ux.dx(1)) * dx + \
             ux * (h + z) * ufl.dot(ufl.as_vector([ux, uy]), n) * v_ux * ds(2) + \
             g * h * h.dx(0) * v_ux * dx + g * h * z.dx(0) * v_ux * dx + \
             g * mannings**2 * (abs(ux) * ux / h_corrected**(1/3)) * v_ux * dx + \
             ux * (h + z) * ufl.dot(ufl.as_vector([ux, uy]), n) * v_ux * ds(1) 
    
    momentum_y = ((h+z)*uy - (h_prev+z)*uy_prev) / dt * v_uy * dx - \
             (uy * (h + z) * ux * v_uy.dx(0) + uy * (h + z) * uy * v_uy.dx(1)) * dx + \
             uy * (h + z) * ufl.dot(ufl.as_vector([ux, uy]), n) * v_uy * ds(2) + \
             g * h * h.dx(1) * v_uy * dx + g * h * z.dx(1) * v_uy * dx + \
             g * mannings**2 * (abs(uy) * uy / h_corrected**(1/3)) * v_uy * dx + \
             uy * (h + z) * ufl.dot(ufl.as_vector([ux, uy]), n) * v_uy * ds(1)  
    return continuity + momentum_x + momentum_y 

# collapsing the spaces for later use in time loop
V0, V0_to_W0 = V.sub(0).collapse()
V1, V1_to_W1 = V.sub(1).collapse()
V2, V2_to_W2 = V.sub(2).collapse()

# boundaries
def open_boundary(x):
    return np.isclose(x[0], 0)
def wall_boundary(x):
    return np.logical_not(np.isclose(x[0], 0)) | np.isclose(x[1], 0) | np.isclose(x[1], Ly)

# boundary facets extraction
boundary_facets_open = locate_entities_boundary(domain, domain.topology.dim - 1, open_boundary)
boundary_facets_wall = locate_entities_boundary(domain, domain.topology.dim - 1, wall_boundary)

# creating tags for corr. facets
facet_tags = mesh.meshtags(domain, domain.topology.dim - 1, 
                      np.concatenate([boundary_facets_wall, boundary_facets_open]),
                      np.concatenate([np.full(len(boundary_facets_wall), 1), 
                                      np.full(len(boundary_facets_open), 2)]))

# for boundary integration setting ds
ds = Measure("ds", domain=domain, subdomain_data=facet_tags)

# for tidal wave at each time step
tidal_value = fem.Constant(domain, ScalarType(0.0))
def update_tidal_value(t):
    tidal_value.value = (
        np.tanh(2.0 * t / (86400.0 * 2.0)) * 2.0 * np.cos(t * (2.0 * np.pi / (12.0 * 60 * 60)) - (90 * np.pi / 180))
    )
update_tidal_value(t=0.0) # initialize the tidal wave

# problem set up
problem = NonlinearProblem(weak_form(W, W_prev, V_test, wd=True), W)
solver = NewtonSolver(MPI.COMM_WORLD, problem)
solver.rtol = 1e-6 
solver.atol = 1e-6  
solver.max_it = 500
solver.error_on_nonconvergence = True  

# set up for plotting the mesh
cells, cell_types, points = vtk_mesh(domain)
pv_mesh = pv.UnstructuredGrid(cells, cell_types, points)
plotter = pv.Plotter()

# time steps and staions to store h, ux, uy values
desired_steps = [i for i in range(0, 9000000)]
time_steps = []
station_coords = [(9000, 3650), (11000, 3650), (13500, 3650)]  
stations = {coord: None for coord in station_coords}  
for coord in station_coords:
    x_coord, y_coord = coord  
    distances = np.sqrt((x[:, 0] - x_coord)**2 + (x[:, 1] - y_coord)**2)  
    station_idx = np.argmin(distances)  
    stations[coord] = station_idx  
print("Nearest points in the mesh for each station:", stations)

# dictionary to store the stations data
time_series_data = {coord: {"time": [], "h": [], "ux": [], "uy": []} for coord in station_coords}

# time loop
for step in range(num_steps+1):
    residual = mass_conservation_residual(h, ux, uy, h_prev, ux_prev, uy_prev, Lx / nx, Ly / ny, dt, maps)
    print(f"Step {step}, Mass Conservation Residual: {np.linalg.norm(residual)}")
    t = step * dt
    time_steps.append(t)
    update_tidal_value(t)   # update tidal wave for each time
    if step in desired_steps:
        # creating individual functions to store current individual values to plot
        z_values_func = fem.Function(V0)
        h_values_func = fem.Function(V0)
        ux_values_func = fem.Function(V1)
        uy_values_func = fem.Function(V2)
        z_values_func.x.array[:] = z.x.array[:] # store bathymetry
        print('z_values_func', z_values_func)
        h_values = h.x.array[maps[0]]  # extract h
        h_values_func.x.array[:] = h_values  # store h 
        ux_values = ux.x.array[maps[1]] # extract ux
        ux_values_func.x.array[:] = ux_values # store ux
        uy_values = uy.x.array[maps[2]] # extract uy
        uy_values_func.x.array[:] = uy_values # store uy
        print('h_values', h_values)
        print('ux', type(ux_values))
        print('z', type(z))
        print('h_values_func', type(h_values_func))
        print('h', type(h))

        # find h, ux, uy values at vertices to plot 
        z_values_at_vertices = find_values_at_vertices(z_values_func, domain)
        z_values_array = np.concatenate(z_values_at_vertices)
        h_values_at_vertices = find_values_at_vertices(h_values_func, domain)
        h_values_array = np.concatenate(h_values_at_vertices)
        ux_values_at_vertices = find_values_at_vertices(ux_values_func, domain)
        ux_values_array = np.concatenate(ux_values_at_vertices)
        uy_values_at_vertices = find_values_at_vertices(uy_values_func, domain)
        uy_values_array = np.concatenate(uy_values_at_vertices)

        # storing the reqd. values into pv_mesh for plotting
        pv_mesh.point_data["Bed Elevation (z)"] = z_values_array
        pv_mesh.point_data["Water Depth (h)"] = h_values_array
        pv_mesh.point_data["Bed Elevation (z) + Water Depth (h)"] = z_values_array + h_values_array
        pv_mesh.point_data["Discharge (qx)"] = h_values_array * ux_values_array
        pv_mesh.point_data["Discharge (qy)"] = h_values_array * uy_values_array
        pv.start_xvfb()

        # warping by bed elevation
        warped = pv_mesh.warp_by_scalar("Bed Elevation (z)")
        plotter.add_mesh(
            warped, 
            scalars="Bed Elevation (z) + Water Depth (h)", # plotting total elevation
            cmap="viridis", 
            show_edges=True, 
            show_scalar_bar=True,
            scalar_bar_args={"title": "Water Depth (h)"}
        )
        arrow_magnitude = 0.5
        # arrows = np.column_stack((ux_values, uy_values, np.zeros_like(ux_values))) * arrow_magnitude
        plotter.add_text(f"Time: {t}s", position="upper_left", font_size=12)
        plotter.show()
        filename = f"frames/frame_{step:04d}.png"
        plotter.screenshot(filename)
        plotter.clear()
    # storing the station data
    for coord, idx in stations.items():
        time_series_data[coord]["time"].append(t)
        time_series_data[coord]["h"].append(h.x.array[maps[0]][idx])
        time_series_data[coord]["ux"].append(ux.x.array[maps[1]][idx])
        time_series_data[coord]["uy"].append(uy.x.array[maps[2]][idx])
    # solve and update
    solver.solve(W)
    W.x.scatter_forward()
    h_prev.x.array[:] = h.x.array[:]
    ux_prev.x.array[:] = ux.x.array[:]
    uy_prev.x.array[:] = uy.x.array[:]