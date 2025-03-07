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

    Saves each station's time series as a PNG image.
    """

    print('Creating station time series plots')

    # Iterate through each station coordinate
    for coord in station_coords:
        plt.figure(figsize=(23, 20))
        if coord == (9000, 3650):
            # Plot surface elevation (h)
            plt.subplot(3, 1, 1)
            plt.plot(time_series_data[coord]["time"], time_series_data[coord]["h"], label="Surface Elevation (h)")
            plt.title(f"Station at x = {coord} m")
            plt.xlabel("Time (s)")
            plt.ylabel("Elevation (m)")
            # plt.yticks(np.arange(np.floor(min(time_series_data[coord]["h"])), 
            #                      np.ceil(max(time_series_data[coord]["h"])), 0.5))  # Set y-ticks at 0.5 intervals
            plt.yticks(np.arange(0, 2.9, 0.5))
            plt.grid(axis='y')  # Add grid lines on y-axis
            plt.legend()
    
            # Plot velocity in x-direction (ux)
            plt.subplot(3, 1, 2)
            plt.plot(time_series_data[coord]["time"], time_series_data[coord]["ux"], label="Velocity (ux)", color="orange")
            plt.xlabel("Time (s)")
            plt.ylabel("Velocity (m/s)")
            # plt.yticks(np.arange(np.floor(min(time_series_data[coord]["ux"])) + 0.4, 
            #                      np.ceil(max(time_series_data[coord]["ux"])) + 0.1, 0.2))  # Set y-ticks at 0.2 intervals
            plt.yticks(np.arange(-0.2, 0.49, 0.2))
            plt.grid(axis='y')
            plt.legend()
    
            # Plot velocity in y-direction (uy)
            plt.subplot(3, 1, 3)
            plt.plot(time_series_data[coord]["time"], time_series_data[coord]["uy"], label="Velocity (uy)", color="green")
            plt.xlabel("Time (s)")
            plt.ylabel("Velocity (m/s)")
            plt.yticks(np.arange(-0.6, 0.51, 0.1))  # Set y-ticks at 0.1 intervals
            plt.grid(axis='y')
            plt.legend()
    
            # Adjust layout and save the plot
            plt.tight_layout()
            plt.savefig(f"station_{coord}_timeseries_exp_without_bathymetry_CG.png")
            plt.show()
        if coord == (11000, 3650):
            # Plot surface elevation (h)
            plt.subplot(3, 1, 1)
            plt.plot(time_series_data[coord]["time"], time_series_data[coord]["h"], label="Surface Elevation (h)")
            plt.title(f"Station at x = {coord} m")
            plt.xlabel("Time (s)")
            plt.ylabel("Elevation (m)")
            # plt.yticks(np.arange(np.floor(min(time_series_data[coord]["h"])), 
            #                      np.ceil(max(time_series_data[coord]["h"])), 0.5))  # Set y-ticks at 0.5 intervals
            plt.yticks(np.arange(-0.5, 3, 0.5))
            plt.grid(axis='y')  # Add grid lines on y-axis
            plt.legend()
    
            # Plot velocity in x-direction (ux)
            plt.subplot(3, 1, 2)
            plt.plot(time_series_data[coord]["time"], time_series_data[coord]["ux"], label="Velocity (ux)", color="orange")
            plt.xlabel("Time (s)")
            plt.ylabel("Velocity (m/s)")
            # plt.yticks(np.arange(np.floor(min(time_series_data[coord]["ux"])) + 0.4, 
            #                      np.ceil(max(time_series_data[coord]["ux"])) + 0.1, 0.2))  # Set y-ticks at 0.2 intervals
            plt.yticks(np.arange(-0.2, 0.29, 0.2))
            plt.grid(axis='y')
            plt.legend()
    
            # Plot velocity in y-direction (uy)
            plt.subplot(3, 1, 3)
            plt.plot(time_series_data[coord]["time"], time_series_data[coord]["uy"], label="Velocity (uy)", color="green")
            plt.xlabel("Time (s)")
            plt.ylabel("Velocity (m/s)")
            plt.yticks(np.arange(-0.4, 0.29, 0.1))  # Set y-ticks at 0.1 intervals
            plt.grid(axis='y')
            plt.legend()
    
            # Adjust layout and save the plot
            plt.tight_layout()
            plt.savefig(f"station_{coord}_timeseries_exp_without_bathymetry_CG.png")
            plt.show()
        if coord == (13500, 3650):
            # Plot surface elevation (h)
            plt.subplot(3, 1, 1)
            plt.plot(time_series_data[coord]["time"], time_series_data[coord]["h"], label="Surface Elevation (h)")
            plt.title(f"Station at x = {coord} m")
            plt.xlabel("Time (s)")
            plt.ylabel("Elevation (m)")
            # plt.yticks(np.arange(np.floor(min(time_series_data[coord]["h"])), 
            #                      np.ceil(max(time_series_data[coord]["h"])), 0.5))  # Set y-ticks at 0.5 intervals
            plt.yticks(np.arange(0, 2.59, 0.5))
            plt.grid(axis='y')  # Add grid lines on y-axis
            plt.legend()
    
            # Plot velocity in x-direction (ux)
            plt.subplot(3, 1, 2)
            plt.plot(time_series_data[coord]["time"], time_series_data[coord]["ux"], label="Velocity (ux)", color="orange")
            plt.xlabel("Time (s)")
            plt.ylabel("Velocity (m/s)")
            # plt.yticks(np.arange(np.floor(min(time_series_data[coord]["ux"])) + 0.4, 
            #                      np.ceil(max(time_series_data[coord]["ux"])) + 0.1, 0.2))  # Set y-ticks at 0.2 intervals
            plt.yticks(np.arange(-0.2, 1.0, 0.1))
            plt.grid(axis='y')
            plt.legend()
    
            # Plot velocity in y-direction (uy)
            plt.subplot(3, 1, 3)
            plt.plot(time_series_data[coord]["time"], time_series_data[coord]["uy"], label="Velocity (uy)", color="green")
            plt.xlabel("Time (s)")
            plt.ylabel("Velocity (m/s)")
            # plt.yticks(np.arange(np.floor(min(time_series_data[coord]["uy"])) + 0.7, 
            #                      np.ceil(max(time_series_data[coord]["uy"])) - 0.8, 0.1))  # Set y-ticks at 0.1 intervals
            plt.yticks(np.arange(-0.1, 0.3, 0.1))
            plt.grid(axis='y')
            plt.legend()
    
            # Adjust layout and save the plot
            plt.tight_layout()
            plt.savefig(f"station_{coord}_timeseries_exp_without_bathymetry_CG.png")
            plt.show()

    print('Successfully created and saved station time series plots')