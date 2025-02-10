import matplotlib.pyplot as plt
import numpy as np

def plot_station_timeseries(time_series_data, station_coords):
    print('Creating stations images')
    # Plot time series for each station
    for coord in station_coords:
        plt.figure(figsize=(23, 20))

        # Plot surface elevation
        plt.subplot(3, 1, 1)
        plt.plot(time_series_data[coord]["time"], time_series_data[coord]["h"], label="Surface Elevation (h)")
        plt.title(f"Station at x = {coord} m")
        plt.xlabel("Time (s)")
        plt.ylabel("Elevation (m)")

        # Set y-ticks at unit intervals and enable grid
        y_min, y_max = plt.ylim()  # Get current y-limits
        plt.yticks(np.arange(np.floor(y_min), np.ceil(y_max) + 0.5, 0.5))  # Set y-ticks at unit intervals
        plt.grid(axis='y')  # Add grid lines for y-axis only
        plt.legend()

        # Plot velocity in x-direction
        plt.subplot(3, 1, 2)
        plt.plot(time_series_data[coord]["time"], time_series_data[coord]["ux"], label="Velocity (ux)", color="orange")
        plt.xlabel("Time (s)")
        plt.ylabel("Velocity (m/s)")

        # Set y-ticks at unit intervals and enable grid
        y_min, y_max = plt.ylim()  # Get current y-limits
        plt.yticks(np.arange(np.floor(y_min) + 0.4, np.ceil(y_max) + 0.1, 0.2))  # Set y-ticks at unit intervals
        plt.grid(axis='y')  # Add grid lines for y-axis only
        plt.legend()

        # Plot velocity in y-direction
        plt.subplot(3, 1, 3)
        plt.plot(time_series_data[coord]["time"], time_series_data[coord]["uy"], label="Velocity (uy)", color="green")
        plt.xlabel("Time (s)")
        plt.ylabel("Velocity (m/s)")

        # Set y-ticks at unit intervals and enable grid
        y_min, y_max = plt.ylim()  # Get current y-limits
        plt.yticks(np.arange(np.floor(y_min) + 0.7, np.ceil(y_max) - 0.8, 0.1))  # Set y-ticks at unit intervals
        plt.grid(axis='y')  # Add grid lines for y-axis only
        plt.legend()

        plt.tight_layout()
        plt.savefig(f"station_{coord}_timeseries_slopeBeach_CG.png")  # Save plot as image
        plt.show()
    print('Successfully created stations images')
