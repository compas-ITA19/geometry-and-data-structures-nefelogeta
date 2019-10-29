def cross_array_of_vectors_py (arr_1, arr_2):
    """Given two arrays of vectors of the same size, perform the cross product elemnt-wise w/o numpy.

    Parameters
    ----------
    arr_1 : array – first list of vectors denoted by their XYZ components.
    arr_2 : array – second list of vectors denoted by their XYZ components.

    Returns
    -------
    list
        XYZ components of the resulting vectors.

    """
    
    import sys

    arr_len = len(arr_1)
    if not arr_len == len(arr_2):
        sys.exit("the two arrays must have the same size.")
    
    x_products=[]
    for e in range(arr_len):
        x_products.append([arr_1[e][1] * arr_2[e][2] - arr_1[e][2] * arr_2[e][1], 
                        arr_1[e][2] * arr_2[e][0] - arr_1[e][0] * arr_2[e][2], 
                        arr_1[e][0] * arr_2[e][1] - arr_1[e][1] * arr_2[e][0]])  
    return x_products


def cross_array_of_vectors_np (arr_1, arr_2):
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

    """Given two arrays of vectors of the same size, perform the cross product elemnt-wise w/ numpy.

    Parameters
    ----------
    arr_1 : array – first list of vectors denoted by their XYZ components.
    arr_2 : array – second list of vectors denoted by their XYZ components.

    Returns
    -------
    list
        XYZ components of the resulting vectors.

    """
    import sys
    import numpy as np

    arr_len = len(arr_1)
    if not arr_len == len(arr_2):
        sys.exit("the two arrays must have the same size.")

    x_products = np.cross(arr_1, arr_2)

    return x_products

def orthonormals_from_two_vectors (u, v):
    """Given two non-parallel vectors, creates a set of three orthonormal vectors.

    Parameters
    ----------
    u : (sequence of float) – XYZ components of the first vector.
    v : (sequence of float) – XYZ components of the second vector.

    Returns
    -------
    list
        The three orthonormal vectors.

    """

    import sys
    import numpy as np
    from compas.geometry import cross_vectors

    k = cross_vectors(u, v)
    # check if vectors are parallel
    if np.linalg.norm(k) == 0.0:
        sys.exit('Input vectors cannot be parallel.')
    k = k/np.linalg.norm(k)
    i = u/np.linalg.norm(u)
    j = cross_vectors(k, i)
    j = j/np.linalg.norm(j)

    return [i, j, k]

def convex_polygon_area(polygon):
    """Compute the area of a convex polygon (on the XY plane).

    Parameters
    ----------
    polygon : sequence
        The XY coordinates of the vertices/corners of the polygon.
        The vertices are assumed to be in order.
        The polygon is assumed to be closed:
        the first and last vertex in the sequence should not be the same.

    Returns
    -------
    float
        The area of the polygon.

    """    

    import math
    import sys
    import numpy as np

    from compas.geometry import cross_vectors_xy

    poly_len = len(polygon)

    if poly_len < 3:
        sys.exit("The polygon is not valid")
    elif poly_len==3:
        v1 = np.subtract(polygon[1], polygon[0])
        v2 = np.subtract(polygon[2], polygon[0])
        poly_area = abs(cross_vectors_xy(v2,v1)[2]/2)
    else:
        x_products = []
        # check if the polygon is convex TODO: revew
        for p in range(poly_len):
            AB = np.subtract(polygon[p-1], polygon[p-2])
            BC = np.subtract(polygon[p], polygon[p-1])
            x_products.append(cross_vectors_xy(AB,BC))
        if not (all(i[2]>0.0 for i in x_products) or all(i[2]<0.0 for i in x_products)):
            sys.exit("the polygon is not convex.")
        # compute internal areas 
        int_areas = []
        for p in range(poly_len-2):
            v1 = np.subtract(polygon[p+1], polygon[0])
            v2 = np.subtract(polygon[p+2], polygon[0])
            int_areas.append(abs(cross_vectors_xy(v2,v1)[2]/2))
        #compute toal area
        poly_area = sum(int_areas)

    return poly_area



if __name__ == "__main__":

    # Test task 1.1
    u = [10, 10.0, 0.0] 
    v = [-10.0, 10.0, 0.0]
    i,j,k = orthonormals_from_two_vectors(u, v)
    print(i,j,k)

    # Test task 1.2
    from compas_plotters import Plotter

    plotter = Plotter()

    my_polygon = [
                [0.0, 0.0],
                [1.0, 0.0],
                [1.0, 0.5],
                [0.5, 2.0],
                [0.0, 1.0]
                ]

    plotter.draw_polygons([{'points': my_polygon}])
    plotter.show()

    print(convex_polygon_area(my_polygon))

    # Test task 1.3
    arr_1 = [[2.0, 0.0, 0.0],[1.0, 0.0, 1.0],[1.0, 0.5, 0.0]]
    arr_2 = [[0.5, 3.0, 0.0],[1.0, 0.0, 0.0],[1.0, 5.0, 2.0]]
    ## pure python
    print(cross_array_of_vectors_py(arr_1, arr_2))
    ## using numpy
    print(cross_array_of_vectors_np(arr_1, arr_2))
