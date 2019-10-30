# from compas.geometry import area_polygon 
# from compas.geometry import centroid_polygon 
# from compas.geometry import centroid_polyhedron 
# from compas.geometry import circle_from_points 
# from compas.geometry import convex_hull 
# from compas.geometry import decompose_matrix 
# from compas.geometry import distance_point_line 
# from compas.geometry import distance_line_line 
# from compas.geometry import intersection_line_line 
# from compas.geometry import intersection_line_triangle 
# from compas.geometry import intersection_plane_plane 
# from compas.geometry import is_coplanar 
# from compas.geometry import is_point_in_triangle 
# from compas.geometry import local_to_world_coords 
# from compas.geometry import matrix_from_basis_vectors 
# from compas.geometry import normal_polygon 
# from compas.geometry import normal_triangle 
# from compas.geometry import offset_line 
# from compas.geometry import orient_points 
# from compas.geometry import orthonormalize_axes 
# from compas.geometry import plane_from_points 
# from compas.geometry import reflect_line_triangle 
# from compas.geometry import volume_polyhedron 

import os 
import compas 
from compas.datastructures import Mesh
from compas_plotters import MeshPlotter
from compas.geometry import cross_vectors
from compas.utilities import pairwise
import numpy as np

def compas_example():
    from random import choice

    import compas

    from compas.utilities import pairwise
    from compas.datastructures import Network
    from compas.topology import dijkstra_path
    from compas_plotters import NetworkPlotter


    # make a network from a sample file

    network = Network.from_obj(compas.get('grid_irregular.obj'))


    # start and end

    leaves = list(network.vertices_where({'vertex_degree': 1}))

    start = end = 0
    while start == end:
        start = choice(leaves)
        end = choice(leaves)

    # construc an adjacency dict
    # add weight to the edges corresponding to their length
    # compute the shortest path

    adjacency = {key: network.vertex_neighbors(key) for key in network.vertices()}

    weight = {(u, v): network.edge_length(u, v) for u, v in network.edges()}
    weight.update({(v, u): weight[(u, v)] for u, v in network.edges()})

    path = dijkstra_path(adjacency, weight, start, end)


    # visualize the result

    plotter = NetworkPlotter(network, figsize=(10, 8), fontsize=6)

    edges = []
    for u, v in pairwise(path):
        if v not in network.edge[u]:
            u, v = v, u
        edges.append([u, v])

    plotter.draw_vertices(
        text={key: key for key in (start, end)},
        facecolor={key: '#ff0000' for key in (path[0], path[-1])},
        radius=0.15
    )

    plotter.draw_edges(
        color={(u, v): '#ff0000' for u, v in edges},
        width={(u, v): 3.0 for u, v in edges},
        text={(u, v): '{:.1f}'.format(weight[(u, v)]) for u, v in network.edges()}
    )

    plotter.show()

def traverse_mesh_straight (mesh, p_start, direction):

    """Traverses a rectangular mesh from boundary to boundary in a "straight" line.

    Parameters
    ----------
    mesh : mesh object
        Mesh on which perform the operation.
    p_start: int
        Index of the vertex on the boundary of the mesh.

    Returns
    -------
    None

    """  

    x_values=[]
    y_values=[]
    for v in mesh.vertices_on_boundary():
        x_values.append(mesh.vertex_coordinates(v)[0])
        y_values.append(mesh.vertex_coordinates(v)[1])
    
    x_min=min(x_values)
    x_max=max(x_values)
    y_min=min(y_values)
    y_max=max(y_values)

    mesh_width = x_max - x_min
    mesh_height = y_max - y_min

    x0_boundary = [key for key in mesh.vertices_where({'x': x_min})]
    x1_boundary = [key for key in mesh.vertices_where({'x': x_max})]
    y0_boundary = [key for key in mesh.vertices_where({'y': y_min})]
    y1_boundary = [key for key in mesh.vertices_where({'y': y_max})]

    mesh_size_x = x1_boundary[0]+1.0
    mesh_sixe_y = (x1_boundary[-1]+1)/(x1_boundary[0]+1.0)

    mesh_vector = [1.0, mesh_size_x, 0.0]
    move_vector = [direction[i]*mesh_vector[i] for i in range(3)] 
    move_step = sum(move_vector)

    if p_start <mesh_size_x:
        x_start_index = int(p_start)
    elif p_start in x0_boundary:
        x_start_index   = 0.0
    elif p_start in x1_boundary:
        x_start_index = 5
    else:
        x_start_index = int(p_start-(mesh_size_x*(mesh_sixe_y-1)))
    
    y_start_index = int((p_start/mesh_size_x))

    iter_vector = [(mesh_size_x-1 - x_start_index), (mesh_sixe_y-1-y_start_index), 0.0]
    max_iter = int(sum([direction[i]*iter_vector[i] for i in range(3)]))

    sequence=[]
    for i in range(max_iter+1):
    
        sequence.append(p_start+move_step*i)

    # visualize the result

    plotter = MeshPlotter(mesh, figsize=(10, 8), fontsize=6)

    edges = []
    for u, v in pairwise(sequence):
        edges.append([u, v])

    plotter.draw_vertices(
        facecolor={key: '#ff0000' for key in (sequence[0], sequence[-1])},
        radius=0.15,
    )

    plotter.draw_edges(
        color={(u, v): '#ff0000' for u, v in edges},
        width={(u, v): 3.0 for u, v in edges},
    )

    plotter.show()




if __name__ == "__main__":

    import random
    
    HERE = os.path.dirname(__file__) 
    DATA = os.path.join(HERE, 'data') 
    FILE = os.path.join(DATA, 'faces.obj')

    my_mesh = Mesh.from_obj(FILE)
    my_p_start = random.choice(my_mesh.vertices_on_boundary())
    directions = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0] ]
    my_dir = random.choice(directions)

    traverse_mesh_straight(my_mesh, my_p_start, my_dir)

    pass

