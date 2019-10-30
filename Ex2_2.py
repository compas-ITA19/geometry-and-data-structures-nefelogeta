import os 
import sys

import compas 
from compas.datastructures import Mesh
from compas_plotters import MeshPlotter
from compas.geometry import cross_vectors
from compas.utilities import pairwise
import numpy as np

def mesh_prop(mesh):
    ''' not used in this exercise'''

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

    return [mesh_size_x, mesh_sixe_y, mesh_height, mesh, x_min, x_max, y_min, y_max, x0_boundary, y0_boundary]


def traverse_mesh_straight (mesh, p_start, direction):

    """Traverses a rectangular mesh from boundary to boundary in a "straight" line.

    Parameters
    ----------
    mesh : mesh object
        Mesh on which perform the operation.
    p_start: int
        Index of the vertex on the boundary of the mesh.
    direction: squence of floats
        X or Y directions along which perform the operation. (positive or negative)

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

    x0_boundary = [key for key in mesh.vertices_where({'x': x_min})]
    x1_boundary = [key for key in mesh.vertices_where({'x': x_max})]

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

    if sum(direction)<0:
        iter_vector = [(x_start_index), (y_start_index), 0.0]
    elif sum(direction)==0:
        sys.exit('direction vector cannot be null')
    else:
        iter_vector = [(mesh_size_x-1 - x_start_index), (mesh_sixe_y-1-y_start_index), 0.0]
    max_iter = int(sum([abs(direction[i]*iter_vector[i]) for i in range(3)]))

    sequence=[]
    for i in range(max_iter + 1):
        sequence.append(p_start+move_step*i)

    # visualize the result

    plotter = MeshPlotter(mesh, figsize=(10, 8), fontsize=6)

    edges = []
    for u, v in pairwise(sequence):
        if sum(direction)>0:
            edges.append([u, v])
        else:
            edges.append([v, u])

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
    directions = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0] ]
    my_dir = random.choice(directions)

    print(my_p_start)
    print(my_dir)
    traverse_mesh_straight(my_mesh, my_p_start, my_dir)

