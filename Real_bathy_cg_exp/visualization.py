import pyvista as pv
from dolfinx.plot import vtk_mesh

def setup_visualization(domain):
    cells, cell_types, points = vtk_mesh(domain)
    pv_mesh = pv.UnstructuredGrid(cells, cell_types, points)
    # Create the plotter
    plotter = pv.Plotter(off_screen=True)  # Enable off-screen rendering
    return pv_mesh, plotter
