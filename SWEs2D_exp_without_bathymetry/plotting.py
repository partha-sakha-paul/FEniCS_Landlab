import matplotlib.pyplot as plt
import numpy as np

def plot_station_timeseries(time_series_data, station_coords):
    """
    Plots time series data for water surface elevation and velocity components
    at specified station locations.

    Parameters:
    - time_series_data: Dictionary containing time series data for each station.
      Expected format: 
      {
        coord1: {"time": [...], "h": [...], "ux": [...], "uy": [...]},
        coord2: {...}, 
        ...
      }
    - station_coords: List of x-coordinates where data is recorded.

    Saves each station's time series as a PNG image with 300 DPI.
    """

    print('Creating station time series plots')

    for coord in station_coords:
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(20, 8))  # Side-by-side plots

        # Surface Elevation (h) on the left
        axes[0].plot(time_series_data[coord]["time"], time_series_data[coord]["h"], label="Surface Elevation (h)")
        axes[0].set_title(f"Station at x = {coord} m - Elevation")
        axes[0].set_xlabel("Time (s)")
        axes[0].set_ylabel("Elevation (m)")
        axes[0].set_yticks(np.arange(0, 3, 0.5))  # Adjust as per data
        axes[0].grid(axis='y')
        axes[0].legend()

        # Velocity components on the right
        axes[1].plot(time_series_data[coord]["time"], time_series_data[coord]["ux"], label="Velocity (ux)", color="orange")
        axes[1].plot(time_series_data[coord]["time"], time_series_data[coord]["uy"], label="Velocity (uy)", color="green")
        axes[1].set_title(f"Station at x = {coord} m - Velocities")
        axes[1].set_xlabel("Time (s)")
        axes[1].set_ylabel("Velocity (m/s)")
        axes[1].grid(axis='y')
        axes[1].legend()

        # Adjust layout and save the plot with 300 DPI
        plt.tight_layout()
        plt.savefig(f"station_{coord}_timeseries_exp_without_bathymetry_CG.png", dpi=300)
        plt.show()

    print('Successfully created and saved station time series plots')
