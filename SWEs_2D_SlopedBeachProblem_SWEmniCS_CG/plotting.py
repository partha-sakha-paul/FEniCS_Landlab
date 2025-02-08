import matplotlib.pyplot as plt
import numpy as np

def plot_station_timeseries(time_series_data, station_coords):
    """
    Plot time series data for multiple stations.

    Parameters:
    - time_series_data (dict): Dictionary containing time-series data for each station.
      Expected format:
      {
          coord_1: {"time": [...], "h": [...], "ux": [...], "uy": [...]},
          coord_2: {...},
          ...
      }
    - station_coords (list): List of station coordinates where data is recorded.

    Generates:
    - Plots of surface elevation (h), velocity in x-direction (ux), and velocity in y-direction (uy).
    - Saves plots as images (e.g., "station_<coord>_timeseries.png").
    """

    print("Creating station time series plots...")

    # Loop over each station coordinate
    for coord in station_coords:
        plt.figure(figsize=(10, 6))

        # Surface Elevation (h)
        plt.subplot(3, 1, 1)
        plt.plot(time_series_data[coord]["time"], time_series_data[coord]["h"], label="Surface Elevation (h)")
        plt.title(f"Station at x = {coord} m")
        plt.xlabel("Time (s)")
        plt.ylabel("Elevation (m)")
        
        # Adjust y-axis ticks for better readability
        y_min, y_max = plt.ylim()
        plt.yticks(np.arange(np.floor(y_min), np.ceil(y_max) + 0.5, 0.5))
        plt.grid(axis="y")
        plt.legend()

        # Velocity in x-direction (ux)
        plt.subplot(3, 1, 2)
        plt.plot(time_series_data[coord]["time"], time_series_data[coord]["ux"], label="Velocity (ux)", color="orange")
        plt.xlabel("Time (s)")
        plt.ylabel("Velocity (m/s)")

        # Adjust y-axis ticks
        y_min, y_max = plt.ylim()
        plt.yticks(np.arange(np.floor(y_min) + 0.4, np.ceil(y_max) + 0.1, 0.2))
        plt.grid(axis="y")
        plt.legend()

        # Velocity in y-direction (uy)
        plt.subplot(3, 1, 3)
        plt.plot(time_series_data[coord]["time"], time_series_data[coord]["uy"], label="Velocity (uy)", color="green")
        plt.xlabel("Time (s)")
        plt.ylabel("Velocity (m/s)")

        # Adjust y-axis ticks
        y_min, y_max = plt.ylim()
        plt.yticks(np.arange(np.floor(y_min) + 0.7, np.ceil(y_max) - 0.8, 0.1))
        plt.grid(axis="y")
        plt.legend()

        # Adjust layout and save figure
        plt.tight_layout()
        plt.savefig(f"station_{coord}_timeseries.png")  # Save plot as an image
        plt.show()

    print("Successfully created station time series plots.")
