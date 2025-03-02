import numpy as np
import matplotlib.pyplot as plt

def load_csv_data():
    """Loads all required CSV data files and returns them as NumPy arrays."""
    FenLand_height_values = np.loadtxt("CG_wse_real_bathy.csv", delimiter=",")  
    FenLand_x_vel_values = np.loadtxt("CG_x_vel_real_bathy.csv", delimiter=",")
    FenLand_y_vel_values = np.loadtxt("CG_y_vel_real_bathy.csv", delimiter=",")

    return (FenLand_height_values, FenLand_x_vel_values, FenLand_y_vel_values)

def generate_time_steps():
    """Generates time steps from 0 to 7 days in intervals of 600s."""
    time_steps = np.arange(0, 7 * 24 * 3600 + 1, 600)  # 0s to 7 days, step 600s
    time_days = time_steps / (24 * 3600)  # Convert seconds to days
    return time_days

def save_plot(time, FenLand_values, ylabel, title, FenLand_color, filename):
    """Plots and saves the graph as an image."""
    plt.figure(figsize=(10, 6))  # Set figure size
    
    plt.plot(time, FenLand_values, label="FEniCSx_Landlab", color=FenLand_color, linestyle="-", linewidth=2, alpha=0.8)

    plt.xlabel("Time (Days)", fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.title(title, fontsize=14)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend(fontsize=10)
    
    plt.savefig(filename, dpi=300, bbox_inches="tight")  # Save image with high resolution
    plt.close()  # Close the plot to free memory

def generate_plots():
    """Loads data, processes it, and saves plots."""
    # Load data
    FenLand_height, FenLand_x_vel, FenLand_y_vel = load_csv_data()
    
    # Generate time steps
    time_days = generate_time_steps()

    # Define station labels and variables
    station_labels = ["Station 1 - (6450m, 2470m)", "Station 2 - (6450m, 6180m)", "Station 3 - (6450m, 11110m)"]
    variables = ["Water Surface Elevation (m)", "X Velocity (m/s)", "Y Velocity (m/s)"]
    # Define variable names without units for filenames
    variable_names = ["water_surface_elevation", "x_velocity", "y_velocity"]  # Names without units

    # Organize data arrays
    FenLand_data = [FenLand_height, FenLand_x_vel, FenLand_y_vel]

    # Define colors for each dataset
    FenLand_colors = ["blue", "green", "red"]

    # Loop through each station and save individual plots
    for i in range(3):  # 3 stations
        for j in range(3):  # 3 variables (height, x_vel, y_vel)
            filename = f"station_{i+1}_{variable_names[j]}.png"  # Clean filename
            save_plot(time_days, FenLand_data[j][:1009, i], variables[j], 
                      f"{station_labels[i]} - {variables[j]}", FenLand_colors[i], filename)

    print("All plots have been saved successfully!")
