import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from dolfinx.mesh import Mesh

def assign_bathymetry_to_mesh(domain: Mesh, csv_path: str) -> np.ndarray:
    """
    Assigns bathymetry (elevation) values from DEM data to the mesh points using nearest-neighbor mapping.

    Parameters:
    domain (dolfinx.mesh.Mesh): The computational mesh.
    csv_path (str): Path to the CSV file containing DEM coordinates and elevation.

    Returns:
    np.ndarray: Array of normalized bathymetry values (scaled from 0 to 10).
    """
    # Extract mesh node coordinates
    mesh_coords = domain.geometry.x[:, :2]  # Extract (x, y) coordinates

    # Load CSV DEM data
    df = pd.read_csv(csv_path)

    # Extract x, y, and elevation (Z) values from DEM
    x_values = df["X"].values  # Longitude
    y_values = df["Y"].values  # Latitude
    dem_values = df["Z"].values  # Elevation values

    # Convert DEM longitude-latitude to meters
    mean_latitude = np.mean(y_values)
    meters_per_degree_lon = 111320 * np.cos(np.radians(mean_latitude))
    meters_per_degree_lat = 111320

    # Convert DEM points to meters
    dem_points_m = np.column_stack([
        (x_values - x_values.min()) * meters_per_degree_lon,
        (y_values - y_values.min()) * meters_per_degree_lat
    ])

    # Create KDTree from DEM points
    tree = cKDTree(dem_points_m)

    # Find the nearest DEM point for each mesh point
    distances, indices = tree.query(mesh_coords)

    # Assign elevation values to the mesh points
    bathymetry = dem_values[indices]

    # Automatically determine min and max values
    bathymetry_min = bathymetry.min()
    bathymetry_max = bathymetry.max()

    # Normalize bathymetry to range [0, 10]
    normalized_bathymetry = (bathymetry - bathymetry_min) / (bathymetry_max - bathymetry_min) * 5

    return normalized_bathymetry
