import pyvista as pv
from dolfinx.plot import vtk_mesh

def setup_visualization(domain):
    """
    Sets up the PyVista visualization for a given finite element mesh.

    Parameters:
    - domain: The finite element mesh (dolfinx.mesh).

    Returns:
    - pv_mesh: A PyVista UnstructuredGrid representing the mesh.
    - plotter: A PyVista Plotter object for rendering.

    This function extracts mesh data from the dolfinx domain and converts it into a 
    format that PyVista can visualize. The plotter is set up for off-screen rendering.
    """
    # Extract mesh data (connectivity, cell types, and node positions)
    cells, cell_types, points = vtk_mesh(domain)
    
    # Create a PyVista UnstructuredGrid for visualization
    pv_mesh = pv.UnstructuredGrid(cells, cell_types, points)
    
    # Create the plotter with off-screen rendering enabled (for batch processing)
    plotter = pv.Plotter(off_screen=True)

    return pv_mesh, plotter
