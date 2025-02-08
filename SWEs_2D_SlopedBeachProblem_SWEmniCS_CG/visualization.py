import pyvista as pv
from dolfinx.plot import vtk_mesh

def setup_visualization(domain):
    """
    Set up the visualization for the computational domain using PyVista.

    Parameters:
    - domain: The computational mesh (dolfinx.mesh.Mesh).

    Returns:
    - pv_mesh: PyVista UnstructuredGrid representation of the domain.
    - plotter: PyVista Plotter instance for rendering.
    """

    # Extract mesh data (connectivity, cell types, and point coordinates) for visualization
    cells, cell_types, points = vtk_mesh(domain)

    # Create a PyVista unstructured grid from the extracted mesh data
    pv_mesh = pv.UnstructuredGrid(cells, cell_types, points)

    # Create a PyVista plotter with off-screen rendering enabled
    plotter = pv.Plotter(off_screen=True)

    return pv_mesh, plotter
