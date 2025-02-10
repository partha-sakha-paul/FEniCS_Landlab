# import necessary libraries and dependencies
import numpy as np
from dolfinx import fem
from petsc4py.PETSc import ScalarType
from ufl import TrialFunction, TestFunction, dx
import ufl
import pyvista as pv
import matplotlib.pyplot as plt
import os

# import necessary functions
from domain import create_domain, create_function_space
from boundary_conditions import setup_boundary_conditions
from solver import setup_solver
from visualization import setup_visualization
from flux import compute_flux
from mass_conservation import mass_conservation_residual
from plotting import plot_station_timeseries

# making folder to store the plots for making the video
os.makedirs("frames", exist_ok=True)

# Domain and Function Space
Lx, Ly = 13800.0, 7200.0  
nx, ny = 12, 6
domain = create_domain(Lx, Ly, nx, ny)
V = create_function_space(domain)

# parameters and constants
g = 9.81  
mannings = 0.02
t, dt, T = 0.0, 600, 7*24*3600.0
num_steps = int(T / dt)

# collapsing subspaces for getting corresponding dofs
num_subs = V.num_sub_spaces
spaces = []
maps = []
for i in range(num_subs):
    space_i, map_i = V.sub(i).collapse()
    spaces.append(space_i)
    maps.append(map_i)

# functions for current time step
W = fem.Function(V)  
h, ux, uy = W.split()

# creating bathymetry by updating z coordinates
V0, _ = V.sub(0).collapse()  
z = fem.Function(V0) 
z.interpolate(lambda x: 5.0 / 13800 * (13800 - x[0]))
x = domain.geometry.x

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

# building wetting-drying condition
wd_alpha = 0.36
wd_alpha_sq = ScalarType(wd_alpha**2)
def wd_correction(h, alpha_sq):
    """Compute correction to water column height for wetting-and-drying."""
    return 0.5 * (ufl.sqrt(h**2 + alpha_sq) - h)

# initialize trial and test functions for the weak form
U = TrialFunction(V)
V_test = TestFunction(V)

# weak form for cg elements
def weak_form(U, U_prev, V_test, wd=True):
    h, ux, uy = U
    h_prev, ux_prev, uy_prev = U_prev
    v_h, v_ux, v_uy = V_test
    # print("Symbolic h before wd:", h)
    if wd:
        h_corrected = h + wd_correction(h, wd_alpha_sq)
    else:
        h_corrected = h
    # print("Symbolic h after wd:", h_corrected)
    continuity = (h - h_prev) / dt * v_h * dx + \
                 (((h+z) * ux).dx(0) * v_h + ((h+z) * uy).dx(1) * v_h) * dx
    momentum_x = ((h+z)*ux - (h_prev+z)*ux_prev) / dt * v_ux * dx + \
                 ((h+z)*ux**2 + 0.5*g*h*(h+2*z)).dx(0) * v_ux * dx + ((h+z)*ux*uy).dx(1) * v_ux * dx - \
                 g * h * z.dx(0) * v_ux * dx + \
                 g * mannings**2 * ((abs(ux) * ux) / h_corrected**(1/3)) * v_ux * dx
    momentum_y = ((h+z)*uy - (h_prev+z)*uy_prev) / dt * v_uy * dx + \
                 ((h+z)*uy**2 + 0.5*g*h*(h+2*z)).dx(1) * v_uy * dx + ((h+z)*ux*uy).dx(0) * v_uy * dx - \
                 g * h * z.dx(1) * v_uy * dx + \
                 g * mannings**2 * ((abs(uy) * uy) / h_corrected**(1/3)) * v_uy * dx
    return continuity + momentum_x + momentum_y

# Boundary Conditions set up for open and wall boundaries
bcs, tidal_value = setup_boundary_conditions(domain, V, Ly)
def update_tidal_value(t):
    tidal_value.value = (
        np.tanh(2.0 * t / (86400.0 * 2.0)) * 2.0 * np.cos(t * (2.0 * np.pi / (12.0 * 60 * 60)) - (90 * np.pi / 180))
    )
update_tidal_value(t=0.0)   # initialize the tidal value

# Solver with problem
solver = setup_solver(weak_form, W, W_prev, V_test, bcs)

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
    t = step * dt  # Compute current time
    time_steps.append(t)  # Store the time value

    # Update the tidal value at the current time step
    update_tidal_value(t=t)  

    # Perform computations at specific time steps
    if step in desired_steps:
        print(f"Time step {step}, Time {t:.2f} s")
        # Compute mass conservation residuals
        residual = mass_conservation_residual(h, ux, uy, h_prev, ux_prev, uy_prev, Lx / nx, Ly / ny, dt, maps)
        print(f"Mass Conservation Residual: {np.linalg.norm(residual)}")

        # Compute inflow and outflow rates
        inflow, outflow = compute_flux(h, ux, uy, domain)

        # Print flux: inflow and outflow details
        print(f"Inflow rate: {inflow:.6f}, Outflow rate: {outflow:.6f}")
        print(f"Net flux: {inflow - outflow:.6f}")

        # Extract field values from finite element solution
        z_values = z.x.array
        h_values = h.x.array[maps[0]]
        ux_values = ux.x.array[maps[1]]
        uy_values = uy.x.array[maps[2]]
        # print(h.x.array)
        # print(h_values)
        # print(ux_values)
        # print(uy_values)

        # Assign values to PyVista mesh for visualization
        pv_mesh.point_data["Bed Elevation (z)"] = z_values
        pv_mesh.point_data["Water Depth (h)"] = h_values
        pv_mesh.point_data["Bed Elevation (z) + Water Depth (h)"] = z_values + h_values
        pv_mesh.point_data["Discharge (qx)"] = h_values * ux_values  # Compute discharge in x-direction
        pv_mesh.point_data["Discharge (qy)"] = h_values * uy_values  # Compute discharge in y-direction

        # Start PyVista's XVFB (virtual framebuffer) mode
        pv.start_xvfb()

        # Warp the mesh by bed elevation for better visualization
        warped = pv_mesh.warp_by_scalar("Bed Elevation (z)")

        # Add the warped mesh to the plot
        plotter.add_mesh(
            warped,
            scalars="Bed Elevation (z) + Water Depth (h)",  # Use water surface elevation as scalar
            cmap="viridis",
            show_edges=True,
            show_scalar_bar=True,
            scalar_bar_args={"title": "Total Surface Elevation (h + z)"}
        )

        # Compute arrow vectors for velocity visualization
        # arrow_magnitude = 0.5  # Scaling factor for arrows
        # arrows = np.column_stack((ux_values, uy_values, np.zeros_like(ux_values))) * arrow_magnitude

        # Add text annotations for inflow, outflow, and net flux
        plotter.add_text(
            f"Inflow: {inflow:.6f}\nOutflow: {outflow:.6f}\nNet flux: {inflow - outflow:.6f}", 
            position="lower_left", font_size=8)

        # Add step information and residual error text
        plotter.add_text(f"Step: {step}, Time: {t:.2f}s\nMass Conservation Residual: {np.linalg.norm(residual)}", font_size=8)

        # Save the current frame as an image
        filename = f"frames/frame_{step:04d}.png"
        plotter.screenshot(filename)

        plotter.clear()  # Clear the plot for the next iteration

    # Store time series data for each station
    for coord, idx in stations.items():
        time_series_data[coord]["time"].append(t)
        time_series_data[coord]["h"].append(h.x.array[maps[0]][idx])  # Store water depth
        time_series_data[coord]["ux"].append(ux.x.array[maps[1]][idx])  # Store velocity in x-direction
        time_series_data[coord]["uy"].append(uy.x.array[maps[2]][idx])  # Store velocity in y-direction

    # Solve the system for the next time step
    solver.solve(W)

    # Update previous time step values for next iteration
    W.x.scatter_forward()
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
