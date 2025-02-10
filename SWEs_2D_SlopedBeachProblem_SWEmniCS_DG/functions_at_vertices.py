import numpy as np

def find_values_at_vertices(function, mesh):
    """
    Extract function values at the mesh vertices.

    Parameters:
    function : dolfinx.fem.Function
        The function whose values need to be evaluated.
    mesh : dolfinx.mesh.Mesh
        The computational mesh.

    Returns:
    values_at_vertices : list
        A list containing function values at each vertex.
    """
    # Ensure connectivity between vertices and higher-dimensional entities
    mesh.topology.create_connectivity(0, mesh.topology.dim)  # Vertex-to-cell connectivity
    mesh.topology.create_connectivity(0, 1)  # Vertex-to-edge connectivity

    # Get connectivity information: which cells are linked to each vertex
    cells = mesh.topology.connectivity(0, mesh.topology.dim)

    # Extract vertex coordinates
    vertices = mesh.geometry.x
    values_at_vertices = []  # List to store function values at vertices

    # Iterate over each vertex in the mesh
    for vertex_index, vertex_coords in enumerate(vertices):
        # Convert vertex coordinates to a 3D array (assuming 2D mesh, z = 0)
        points = np.array([[vertex_coords[0], vertex_coords[1], 0.0]])

        # Find cells associated with this vertex
        associated_cells = cells.links(vertex_index)

        # Evaluate function value at the vertex using the first associated cell
        if len(associated_cells) > 0:
            value = function.eval(points, [associated_cells[0]])
            values_at_vertices.append(value)

    return values_at_vertices
