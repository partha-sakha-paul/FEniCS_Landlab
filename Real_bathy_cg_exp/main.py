# import necessary libraries and dependencies
from mpi4py import MPI
from petsc4py import PETSc
import numpy as np
from dolfinx import fem, geometry
from petsc4py.PETSc import ScalarType
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from ufl import TestFunction, dx, as_vector
from dolfinx.fem.petsc import assemble_vector
from dolfinx.fem import form
import pyvista as pv
import os

# import necessary functions
from domain import create_domain
from flux import compute_flux
from mass_conservation import mass_conservation_residual
from domain import create_function_space
from functions_at_vertices import find_values_at_vertices
from weak_form import weak_form
from boundary_conditions import initialize_boundary_conditions
from visualization import setup_visualization
import csv

# making folder to store the plots for making the video
os.makedirs("frames", exist_ok=True)

# File names to store the solution variables
wse_file = "CG_wse_real_bathy.csv"
x_vel_file = "CG_x_vel_real_bathy.csv"
y_vel_file = "CG_y_vel_real_bathy.csv"

# Function to append a row to a CSV file
def append_to_csv(file_name, data):
    with open(file_name, mode='a', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(data)  # Write the row of data

# for dem data usage
from extract_dem import extract_dem_to_csv
from convert_dem_to_meters import convert_dem_to_meters
from dolfinx import mesh
from assign_bathymetry import assign_bathymetry_to_mesh

# Define paths
dem_file = "dem_swe.tif"  # add required file
csv_file = "dem_coordinates_for_swes.csv"

# Extract DEM data to CSV
extract_dem_to_csv(dem_file, csv_file)

# Convert DEM coordinates to meters
Lx, Ly = convert_dem_to_meters(csv_file)

# Domain and Function Space
nx, ny = 6, 12
domain = create_domain(Lx, Ly, nx, ny)
V = create_function_space(domain)

# Assign bathymetry values to the mesh
bathymetry = assign_bathymetry_to_mesh(domain, csv_file)

# Now you can use bathymetry in further computations
print("Assigned bathymetry values:", bathymetry.shape)

# parameters and constants
g = 9.81  
mannings = 0.02
t, dt, T = 0.0, 600, 600*7*24*6
num_steps = int(T / dt)

# collapsing subspaces for getting corresponding dofs
num_subs = V.num_sub_spaces
spaces = []
maps = []
for i in range(num_subs):
    space_i, map_i = V.sub(i).collapse()
    spaces.append(space_i)
    maps.append(map_i)
V0, V1 = spaces[0], spaces[1]

# Creating bathymetry by updating z coordinates
z = fem.Function(V0) 
# z.interpolate(lambda x: 5.0 / 13800 * (13800 - x[0]))
z.x.array[:] = bathymetry

# Getting the coordinates of the mesh to get the station coordinates
x = domain.geometry.x

# functions for current time step
W = fem.Function(V)
H, uxuy = W.split()

# initializing current functions with zeros as dry bed and still water
H.interpolate(
    fem.Expression(z, V.sub(0).element.interpolation_points)
)
uxuy.interpolate(
    fem.Expression(
        as_vector([fem.Constant(domain, ScalarType(0.0)),fem.Constant(domain, ScalarType(0.0))]), V.sub(1).element.interpolation_points
    )
)
# print(uxuy.x.array[:][maps[1]])

# functions for previous time step
W_prev = fem.Function(V)
H_prev, uxuy_prev = W_prev.split()

# initializing previous functions with zeros as dry bed and still water
H_prev.interpolate(
    fem.Expression(z, V.sub(0).element.interpolation_points)
)
uxuy_prev.interpolate(
    fem.Expression(
        as_vector([fem.Constant(domain, ScalarType(0.0)),fem.Constant(domain, ScalarType(0.0))]), V.sub(1).element.interpolation_points
    )
)

# functions for previous of previous time step: required for BDF2
W_prevprev = fem.Function(V)  
H_prevprev, uxuy_prevprev = W_prevprev.split()

# functions for open boundary conditions and initialization as of dry bed and still water
W_open = fem.Function(V)
H_open, uxuy_open = W_open.split()
# H_open.interpolate(
#     fem.Expression(z, V.sub(0).element.interpolation_points)
# )
W_open.sub(0).interpolate(
    fem.Expression(z, V.sub(0).element.interpolation_points)
)
W_open.sub(1).interpolate(
    lambda x: np.vstack([np.zeros(x.shape[1]),np.zeros(x.shape[1])])
)

# measure ds for boundary integration using facet tags and extracting dofs of h in open boundary for tidal updating
ds, sub_dofs_h_open_left = initialize_boundary_conditions(V, W_open, domain, Lx, Ly)
# print(sub_dofs_h_open_left)
# Testfunction for the Weak form
V_test = TestFunction(V)

# Visualization set up
pv_mesh, plotter = setup_visualization(domain)

# Desired steps to plot
desired_steps = [i for i in range(0, 9000000)]

# Define station coordinates for extracting time series data
stations = np.array([
    [6450, 2470, 0], [6450, 6180, 0], [6450, 11110, 0]
])

# Create a bounding box tree for efficient spatial queries in the given domain.
bb_tree = geometry.bb_tree(domain, domain.topology.dim)

# Lists to store cell indices and corresponding points that are on the current process.
cells = []
points_on_proc = []

# Compute candidate cells that may contain the given station points.
cell_candidates = geometry.compute_collisions_points(bb_tree, stations)

# Identify the exact cells that contain the station points.
colliding_cells = geometry.compute_colliding_cells(domain, cell_candidates, stations)

# List to store indices of stations that are found within the domain.
station_index = []

# Iterate through each station point to check if it belongs to any cell.
for i, point in enumerate(stations):
    if len(colliding_cells.links(i)) > 0:  # If there is at least one colliding cell
        points_on_proc.append(point)  # Store the station point
        cells.append(colliding_cells.links(i)[0])  # Store the first colliding cell
        station_index.append(i)  # Store the index of the station

# Evaluate the bathymetry (z-values) at the identified station points.
station_bathy = z.eval(points_on_proc, cells)

# Convert the list of valid station points into a NumPy array for efficient numerical operations.
points_on_proc = np.array(points_on_proc, dtype=np.float64)

# Initialize global min/max before the simulation loop
global_H_min, global_H_max = float("inf"), float("-inf")

# Iterate over time steps
for step in range(num_steps + 1):
    t = (step+1) * dt
    # print(W_open.x.array)

    # function for open bcs: to extract the values of bed at open boundary
    z_open = W_open.sub(0)
    z_open.interpolate(
        fem.Expression(z, V.sub(0).element.interpolation_points)
    )
    # print('z.x.array[sub_dofs_h_open_left]', z_open.x.array[sub_dofs_h_open_left])
    # print(W_open.x.array)

    # Update open boundary condition at left boundary
    H_open.x.array[sub_dofs_h_open_left] = z_open.x.array[sub_dofs_h_open_left] + (
        np.tanh(2.0 * t / (86400.0 * 2.0)) * 2.0 * np.cos(t * (2.0 * np.pi / (12.0 * 60 * 60)) - (90 * np.pi / 180))
    )
    W.x.array[sub_dofs_h_open_left] = W_open.x.array[sub_dofs_h_open_left]
    # print(W_open.x.array)

    # Set theta1 parameter based on step index
    theta1 = 0.0 if step < 2 else 1.0

    # Compute weak form
    weak_result = weak_form(W, W_prev, W_prevprev, V_test, W_open, domain, z, dt, ds, g, mannings, step, wd_alpha=0.36, theta1=theta1, wd=True)
    print('weak norm1', np.linalg.norm(assemble_vector(form(weak_result)).array))

    # Define the nonlinear problem
    problem = NonlinearProblem(weak_result, W)

    # Setup Newton solver
    solver = NewtonSolver(MPI.COMM_WORLD, problem)
    solver.rtol = 1e-6
    solver.atol = 1e-6
    solver.max_it = 500
    solver.error_on_nonconvergence = True
    solver.report = True

    # Assemble and check the residual vector
    F_form = form(weak_result)
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
    
    # Compute mass conservation residuals
    residual = mass_conservation_residual(H, uxuy, H_prev, uxuy_prev, Lx / nx, Ly / ny, dt, maps)
    
    # print(f"Residual Norm: {residual_norm}")
    if step in desired_steps:
        # Compute inflow and outflow fluxes
        inflow, outflow = compute_flux(H, uxuy, domain)
        print(f"Time step {step}, Time {t:.2f} s")
        print(f"Mass Conservation Residual: {np.linalg.norm(residual)}")
        print(f"Residual Norm with Newton Solver: {residual_norm}")
        print(f"Inflow rate: {inflow:.6f}, Outflow rate: {outflow:.6f}")
        print(f"Net flux: {inflow - outflow:.6f}")
        
        # Create FEM function objects to store solution values
        z_values_func = fem.Function(V0)
        H_values_func = fem.Function(V0)
    
        # Copy solution values into the function objects
        z_values_func.x.array[:] = z.x.array[:]
        H_values_func.x.array[:] = H.x.array[maps[0]]
    
        # Extract values at vertices
        z_values_at_vertices = find_values_at_vertices(z_values_func, domain)
        z_values_array = np.concatenate(z_values_at_vertices)
    
        H_values_at_vertices = find_values_at_vertices(H_values_func, domain)
        H_values_array = np.concatenate(H_values_at_vertices)
        
        H_values = W.sub(0).eval(points_on_proc, cells)
        h_values = H_values - station_bathy
        u_values = W.sub(1).eval(points_on_proc, cells)
        u_values = np.hstack([h_values,u_values])
        print('u_values',u_values)
        
        # Extract columns for storing the station data in csv files
        height_values = u_values[:, 0]  # First column (Height)
        x_vel_values = u_values[:, 1]   # Second column (X velocity)
        y_vel_values = u_values[:, 2]   # Third column (Y velocity)
    
        # Append data to CSV files
        append_to_csv(wse_file, height_values)
        append_to_csv(x_vel_file, x_vel_values)
        append_to_csv(y_vel_file, y_vel_values)
        
        # Update global min/max with current time step data
        global_H_min = min(global_H_min, np.min(H_values_array))
        global_H_max = max(global_H_max, np.max(H_values_array))
    
        # Assign computed values to the PyVista mesh
        pv_mesh.point_data["Bed Elevation (z)"] = z_values_array
        pv_mesh.point_data["Total Surface Elevation (H)"] = H_values_array
        # pv_mesh.point_data["Discharge (qx)"] = H_values_array * ux_values_array
        # pv_mesh.point_data["Discharge (qy)"] = H_values_array * uy_values_array
        
        # Start virtual framebuffer for rendering (useful for headless environments)
        pv.start_xvfb()
        
        # Define fixed color limits
        vmin, vmax = 0, 5.5  # Set min and max values to manage the colorbar

        # Warp mesh by bed elevation for visualization
        warped = pv_mesh.warp_by_scalar("Bed Elevation (z)")
        
        # Add mesh to plotter with water depth visualization
        plotter.add_mesh(
            warped, 
            scalars="Total Surface Elevation (H)", 
            cmap="viridis", 
            show_edges=True, 
            show_scalar_bar=True,
            scalar_bar_args={"title": "Total Surface Elevation (H = h+z): CG"},
            clim=[vmin, vmax]  # Set fixed range
        )
        
        # Add text annotations for inflow, outflow, and net flux
        plotter.add_text(
            f"Inflow: {inflow:.6f}\nOutflow: {outflow:.6f}\nNet flux: {inflow - outflow:.6f}", 
            position="lower_left", font_size=8)

        # Add step information and residual error text
        plotter.add_text(f"Step: {step}, Time: {t:.2f}s\nMass Conservation Residual: {np.linalg.norm(residual)}", font_size=8)

        # Save screenshot of the current frame
        filename = f"frames/frame_{step:04d}.png"
        plotter.screenshot(filename)
        
        # Clear the plotter for the next frame
        plotter.clear()
    
    # Solve the nonlinear problem
    solver.solve(W)
    W.x.scatter_forward()
    
    # Update previous time-step values
    W_prevprev.x.array[:] = W_prev.x.array[:]
    W_prev.x.array[:] = W.x.array[:]
    
# Print or store the values for later analysis
print(f"Overall Min h = {global_H_min}, Overall Max h = {global_H_max}")

# Call the function to plot time series data for each station
from plotting import generate_plots

# Run the plot generation process
generate_plots()

# Import video generation script
import generate_video

# Call the function to generate a video from the saved frames
generate_video.generate_video()

print("Simulation complete. Images and Video created succesfully!!! Please check the current folder to see the results. Thank you!!")