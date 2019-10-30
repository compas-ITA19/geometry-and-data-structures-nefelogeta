import os 
import sys

import compas 
from compas.datastructures import Mesh
from compas_plotters import MeshPlotter
from compas.geometry import cross_vectors
from compas.utilities import pairwise
from compas.topology import dijkstra_path
import numpy as np
import random


def traverse_mesh(mesh, p_start):
    """Traverses a quad mesh from boundary to boundary in a "straight" line.
    If a corner point is selected as starting point a random choice is performed about
    the direction of the path.

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

    sequence = [p_start]
    neighbors = mesh.vertex_neighbors(p_start)

    current = p_start
    if len(neighbors)==2:
        p = random.choice(neighbors)
        print(neighbors, p)
        previous, current = current, p
        while True:
            sequence.append(current)
            if mesh.vertex_degree(current)==2:
                break
            neighbors = mesh.vertex_neighbors(current, ordered=True)
            previous = current
            for n in neighbors:
                if (mesh.is_vertex_on_boundary(n) and  n != sequence[-2]):
                    current = n
    else:
        for p in neighbors:
            if not mesh.is_vertex_on_boundary(p):
                previous, current = current, p
                break
        while True:
            sequence.append(current)
            if mesh.is_vertex_on_boundary(current):
                break
            neighbors = mesh.vertex_neighbors(current, ordered=True)
            i = neighbors.index(previous)
            previous, current = current, neighbors[i - 2]

    print(sequence)

    edges = []
    plotter = MeshPlotter(mesh, figsize=(10, 8), fontsize=6)

    for u, v in pairwise(sequence):
        if u<v:
            edges.append([u, v])
        else:
            edges.append([v, u])    
    plotter.draw_vertices(
                        facecolor={key: '#ff0000' for key in sequence},
                        radius=0.025,
                        )

    plotter.draw_edges(
                    color={(u, v): '#ff0000' for u, v in edges},
                    width={(u, v): 4.0 for u, v in edges},
                    )

    plotter.show()

if __name__ == "__main__":

    import random
    
    HERE = os.path.dirname(__file__) 
    DATA = os.path.join(HERE, 'data') 
    # FILE = os.path.join(DATA, 'faces.obj')
    # my_mesh = Mesh.from_obj(FILE)

    FILE = os.path.join(DATA, 'fofin.json')
    my_mesh = Mesh.from_json(FILE)

    my_p_start = random.choice(my_mesh.vertices_on_boundary())

    traverse_mesh(my_mesh, my_p_start)


