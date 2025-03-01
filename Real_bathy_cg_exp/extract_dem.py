import rasterio
import numpy as np
import pandas as pd

def extract_dem_to_csv(dem_path: str, output_csv: str) -> None:
    """
    Extracts DEM data from a raster file and saves it as a CSV file.

    Parameters:
    dem_path (str): Path to the input DEM file.
    output_csv (str): Path to save the extracted CSV file.

    Returns:
    None
    """
    with rasterio.open(dem_path) as src:
        elevation = src.read(1)  # Read elevation values
        transform = src.transform  # Get affine transformation

    # Get dimensions
    rows, cols = elevation.shape

    # Create lists to store coordinates
    x_coords, y_coords, z_coords = [], [], []

    # Extract coordinates
    for i in range(rows):
        for j in range(cols):
            x, y = transform * (j, i)  # Convert pixel indices to coordinates
            z = elevation[i, j]
            if z > -100:  # Filter out invalid values
                x_coords.append(x)
                y_coords.append(y)
                z_coords.append(z)

    # Create a DataFrame
    df = pd.DataFrame({"X": x_coords, "Y": y_coords, "Z": z_coords})
    
    # Save as CSV
    df.to_csv(output_csv, index=False)
    print(f"Extracted DEM coordinates saved to {output_csv}")
