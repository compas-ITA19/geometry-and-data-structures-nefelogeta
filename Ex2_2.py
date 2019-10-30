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

# import os 
# import compas 
# from compas.datastructures import Mesh
# from compas_plotters import MeshPlotter

# HERE = os.path.dirname(__file__) 
# DATA = os.path.join(HERE, 'data') 
# FILE = os.path.join(DATA, 'faces.obj')

# mesh = Mesh.from_obj(FILE)
# plotter = MeshPlotter(mesh, figsize=(8, 5))

# plotter.draw_vertices(facecolor={key: (255, 0, 0) for key in mesh.vertices_on_boundary()},    
#                     text={key: str(mesh.vertex_degree(key)) for key in mesh.vertices()},    
#                     radius=0.2)

# plotter.draw_edges(keys=list(mesh.edges_on_boundary()),  
#                             color=(255, 0, 0))

# plotter.draw_faces(facecolor={key: (150, 255, 150) for key in mesh.faces() if not mesh.is_face_on_boundary(key)})
# plotter.show()

