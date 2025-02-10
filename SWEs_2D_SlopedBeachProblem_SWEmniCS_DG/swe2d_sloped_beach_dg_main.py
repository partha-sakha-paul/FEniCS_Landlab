# import necessary libraries and dependencies
from mpi4py import MPI
from petsc4py import PETSc
import numpy as np
from dolfinx import fem
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from ufl import TestFunction, dx
from dolfinx.fem.petsc import assemble_vector
from dolfinx.fem import form
import pyvista as pv
import os

# import necessary functions
from domain import create_domain
from flux import compute_flux
# from mass_conservation import mass_conservation_residual
from domain import create_function_space
from functions_at_vertices import find_values_at_vertices
from weak_form import weak_form
from boundary_conditions import initialize_boundary_conditions
from visualization import setup_visualization
from plotting import plot_station_timeseries


# making folder to store the plots for making the video
os.makedirs("frames", exist_ok=True)

# parameters and constants
Lx, Ly = 13800.0, 7200.0  
nx, ny = 12, 6
g = 9.81  
mannings = 0.02
t, dt, T = 0.0, 600, 7*24*3600
num_steps = int(T / dt)

# Domain and Function Space
domain = create_domain(Lx, Ly, nx, ny)
V = create_function_space(domain)

# collapsing subspaces for getting corresponding dofs
num_subs = V.num_sub_spaces
spaces = []
maps = []
for i in range(num_subs):
    space_i, map_i = V.sub(i).collapse()
    spaces.append(space_i)
    maps.append(map_i)
V0, V1, V2 = spaces[0], spaces[1], spaces[2]

# Creating bathymetry by updating z coordinates
z = fem.Function(V0) 
z.interpolate(lambda x: 5.0 / 13800 * (13800 - x[0]))

# Getting the coordinates of the mesh to get the station coordinates
x = domain.geometry.x

# functions for current time step
W = fem.Function(V)
h, ux, uy = W.split()

# initializing current functions with zeros as dry bed
h.interpolate(lambda x: np.full_like(x[0], 0.0))
ux.interpolate(lambda x: np.full_like(x[0], 0.0))
uy.interpolate(lambda x: np.full_like(x[0], 0.0))

# functions for previous time step
W_prev = fem.Function(V)
h_prev, ux_prev, uy_prev = W_prev.split()

# initializing previous functions with zeros as dry bed
h_prev.interpolate(lambda x: np.full_like(x[0], 0.0))
ux_prev.interpolate(lambda x: np.full_like(x[0], 0.0))
uy_prev.interpolate(lambda x: np.full_like(x[0], 0.0))

W_prevprev = fem.Function(V)  
h_prevprev, ux_prevprev, uy_prevprev = W_prevprev.split()

W_open = fem.Function(V)
h_open, ux_open, uy_open = W_open.split()
h_open.interpolate(lambda x: np.full_like(x[0], 0.000))
ux_open.interpolate(lambda x: np.full_like(x[0], 0.0))
uy_open.interpolate(lambda x: np.full_like(x[0], 0.0))

ds, sub_dofs_h_open_left = initialize_boundary_conditions(V, W_open, domain)

V_test = TestFunction(V)

# Visualization set up
pv_mesh, plotter = setup_visualization(domain)

# Desired steps to plot
desired_steps = [i for i in range(0, 9000000)]

# Store time steps to plot h, ux, uy graphs
time_steps = []

# Define station coordinates for extracting time series data
station_coords = [(9000, 3650), (11000, 3650), (13500, 3650)]

# Initialize a dictionary to store the nearest mesh index for each station
stations = {coord: None for coord in station_coords}

# Find the nearest mesh index for each station coordinate
for coord in station_coords:
    x_coord, y_coord = coord  # Extract x and y coordinates of the station

    # Compute Euclidean distance from each mesh point to the station
    distances = np.sqrt((x[:, 0] - x_coord)**2 + (x[:, 1] - y_coord)**2)  

    # Find the index of the closest mesh point
    station_idx = np.argmin(distances)  
    stations[coord] = station_idx  # Store the index in the dictionary

# Initialize a dictionary to store time series data for each station
time_series_data = {coord: {"time": [], "h": [], "ux": [], "uy": []} for coord in station_coords}

# Iterate over time steps
for step in range(num_steps + 1):
    t = step * dt
    time_steps.append(t)

    # Update open boundary condition at left boundary
    h_open.x.array[sub_dofs_h_open_left] = (
        np.tanh(2.0 * t / (86400.0 * 2.0)) * 2.0 * np.cos(t * (2.0 * np.pi / (12.0 * 60 * 60)) - (90 * np.pi / 180))
    )

    # Set theta1 parameter based on step index
    theta1 = 0.0 if step < 2 else 1.0

    # Define the nonlinear problem
    problem = NonlinearProblem(
        weak_form(W, W_prev, W_prevprev, V_test, W_open, domain, z, dt, ds, g, mannings, step, wd_alpha=0.36, theta1=theta1, wd=True), W
    )

    # Setup Newton solver
    solver = NewtonSolver(MPI.COMM_WORLD, problem)
    solver.rtol = 1e-6
    solver.atol = 1e-6
    solver.max_it = 500
    solver.error_on_nonconvergence = True
    solver.report = True

    # Assemble and check the residual vector
    F_form = form(
        weak_form(W, W_prev, W_prevprev, V_test, W_open, domain, z, dt, ds, g, mannings, step, wd_alpha=0.36, theta1=theta1, wd=True)
    )
    F_vec = assemble_vector(F_form)

    # Compute and print the residual norm
    residual_norm = np.linalg.norm(F_vec.array)

    # # Optional debugging checks for NaN or Inf values in residual vector
    # if np.isnan(F_vec.array).any():
    #     print("Warning: NaN detected in residual vector!")

    # if np.isinf(F_vec.array).any():
    #     print("Warning: Inf detected in residual vector!")

    # nan_indices = np.where(np.isnan(F_vec.array))[0]
    # if len(nan_indices) > 0:
    #     print(f"Indices with NaN: {nan_indices}")
    #     print(f"Corresponding Residual Values: {F_vec.array[nan_indices]}")
    
    if step in desired_steps:
        # Compute inflow and outflow fluxes
        inflow, outflow = compute_flux(h, ux, uy, domain)
        print(f"Time step {step}, Time {t:.2f} s")
        print(f"Residual Norm: {residual_norm}")
        print(f"Inflow rate: {inflow:.6f}, Outflow rate: {outflow:.6f}")
        print(f"Net flux: {inflow - outflow:.6f}")
    
        # Create FEM function objects to store solution values
        z_values_func = fem.Function(V0)
        h_values_func = fem.Function(V0)
        ux_values_func = fem.Function(V1)
        uy_values_func = fem.Function(V2)
    
        # Copy solution values into the function objects
        z_values_func.x.array[:] = z.x.array[:]
        h_values_func.x.array[:] = h.x.array[maps[0]]
        ux_values_func.x.array[:] = ux.x.array[maps[1]]
        uy_values_func.x.array[:] = uy.x.array[maps[2]]
    
        # Extract values at vertices
        z_values_at_vertices = find_values_at_vertices(z_values_func, domain)
        z_values_array = np.concatenate(z_values_at_vertices)
    
        h_values_at_vertices = find_values_at_vertices(h_values_func, domain)
        h_values_array = np.concatenate(h_values_at_vertices)
        # print('h_values_array',h_values_array)
    
        ux_values_at_vertices = find_values_at_vertices(ux_values_func, domain)
        ux_values_array = np.concatenate(ux_values_at_vertices)
        # print('ux_values_array',ux_values_array)
    
        uy_values_at_vertices = find_values_at_vertices(uy_values_func, domain)
        uy_values_array = np.concatenate(uy_values_at_vertices)
        # print('uy_values_array',uy_values_array)

        # Assign computed values to the PyVista mesh
        pv_mesh.point_data["Bed Elevation (z)"] = z_values_array
        pv_mesh.point_data["Water Depth (h)"] = h_values_array
        pv_mesh.point_data["Bed Elevation (z) + Water Depth (h)"] = z_values_array + h_values_array
        pv_mesh.point_data["Discharge (qx)"] = h_values_array * ux_values_array
        pv_mesh.point_data["Discharge (qy)"] = h_values_array * uy_values_array
        
        # Start virtual framebuffer for rendering (useful for headless environments)
        pv.start_xvfb()
        
        # Warp mesh by bed elevation for visualization
        warped = pv_mesh.warp_by_scalar("Bed Elevation (z)")
        
        # Add mesh to plotter with water depth visualization
        plotter.add_mesh(
            warped, 
            scalars="Bed Elevation (z) + Water Depth (h)", 
            cmap="viridis", 
            show_edges=True, 
            show_scalar_bar=True,
            scalar_bar_args={"title": "Water Depth (h)"}
        )
        
        # Add time text annotation
        plotter.add_text(f"Time: {t}s", position="upper_left", font_size=12)
        
        # Save screenshot of the current frame
        filename = f"frames/frame_{step:04d}.png"
        plotter.screenshot(filename)
        
        # Clear the plotter for the next frame
        plotter.clear()
        
    # Store time-series data for selected stations
    for coord, idx in stations.items():
        time_series_data[coord]["time"].append(t)
        time_series_data[coord]["h"].append(h.x.array[maps[0]][idx])
        time_series_data[coord]["ux"].append(ux.x.array[maps[1]][idx])
        time_series_data[coord]["uy"].append(uy.x.array[maps[2]][idx])
    
    # Solve the nonlinear problem
    solver.solve(W)
    W.x.scatter_forward()
    
    # Update previous time-step values
    h_prevprev.x.array[:] = h_prev.x.array[:]
    ux_prevprev.x.array[:] = ux_prev.x.array[:]
    uy_prevprev.x.array[:] = uy_prev.x.array[:]
    
    h_prev.x.array[:] = h.x.array[:]
    ux_prev.x.array[:] = ux.x.array[:]
    uy_prev.x.array[:] = uy.x.array[:]

# Call the function to plot time series data for each station
plot_station_timeseries(time_series_data, station_coords)

# Import video generation script
import generate_video

# Call the function to generate a video from the saved frames
generate_video.generate_video()

print("Simulation complete. Images and Video created succesfully!!! Please check the current folder to see the results. Thank you!!")
