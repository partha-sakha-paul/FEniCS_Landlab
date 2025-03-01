import numpy as np
import pandas as pd

def convert_dem_to_meters(csv_path: str):
    """
    Converts DEM longitude-latitude coordinates to meters.

    Parameters:
    csv_path (str): Path to the input CSV file containing DEM coordinates.

    Returns:
    tuple: (Lx, Ly) domain size in meters
    """
    df = pd.read_csv(csv_path)

    # Extract x, y coordinates
    x_values = df["X"].values  # Longitude
    y_values = df["Y"].values  # Latitude

    # Find the domain size in degrees
    Lx_deg = x_values.max() - x_values.min()  # Longitude difference
    Ly_deg = y_values.max() - y_values.min()  # Latitude difference

    # Convert degrees to meters
    mean_latitude = np.mean(y_values)  # Use the mean latitude for conversion
    Lx = Lx_deg * (111320 * np.cos(np.radians(mean_latitude)))  # Convert longitude to meters
    Ly = Ly_deg * 111320  # Convert latitude to meters

    print(f"Domain size: Lx = {Lx:.2f} meters, Ly = {Ly:.2f} meters")
    return Lx, Ly
