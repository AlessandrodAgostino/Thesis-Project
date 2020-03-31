import numpy as np
import itertools
import matplotlib.pyplot as plt
import time
import pandas as pd

from scipy.spatial import Voronoi, ConvexHull, cKDTree
from SALib.sample import saltelli

from D3branch import *

def _inside_boundaries(vor, reg, boundaries):
    """
    This function check if a 'region' of a given 'Voronoi' tassellation lies
    inside some given 'boundaries'.

    Parameters:
    vor    = scipy.spatial.Voronoi object
    reg    = subset of vor.vertices to check
    boundaries = [ [min_x, max_x], [min_y, max_y], [min_z, max_z] ]

    Return: True if 'all' the verticies lies inside the region. False otherwise.
    """
    if -1 in reg: return False
    if len(reg)>1:
        reg = np.asarray([vor.vertices[pt] for pt in reg])
        if all(reg[:,0] > boundaries[0][0]) and all(reg[:,0] < boundaries[0][1]):
            if all(reg[:,1] > boundaries[1][0]) and all(reg[:,1] < boundaries[1][1]):
                if all(reg[:,2] > boundaries[2][0]) and all(reg[:,2] < boundaries[2][1]):
                    return True
        return False
    else: return False

def _plane_y_intersection(p1, p2, k=0):
    """
    Return the intersection between the plane of equation y = k and the line
    crossing 'p1' and 'p2'.
    """
    x = p1[0] + (k - p1[1]) / (p2[1] - p1[1]) * (p2[0] - p1[0])
    z = p1[2] + (k - p1[1]) / (p2[1] - p1[1]) * (p2[2] - p1[2])
    return np.asarray((x,k,z))

def _box_and_spheres(ramification):
    """
    Given a ramification this function returns the list of free end's spheres,
    the radius of those spheres and the box that contains them all with a narrow
    padding, proportional to the spheres' radius.
    """
    #Extracting free end's spheres and radius
    spheres  = [] #List of spheres
    sph_rad  = 0
    max_iter = np.log2((len(ramification)+1)) - 1
    l = ramification[0].l

    for br in ramification:
        if br.iter_lev == max_iter:
            cent = br.pos + vector(*(br.length * br.drct))
            spheres.append(np.asarray((cent.x, cent.y, cent.z)))
            if not sph_rad: sph_rad = br.length

    #Interesting spheres. Those around x,y-plane
    int_spheres = [sph for sph in spheres if (sph[1] > l - 2*sph_rad) & (sph[1] < l + 2*sph_rad)]

    #Boundaries for random sampling in a volume with a padding proportional to spheres' radius
    max_box = np.max(np.asarray(spheres), axis = 0) + sph_rad*2
    min_box = np.min(np.asarray(spheres), axis = 0) - sph_rad*2

    boundaries = [ [min_box[0], max_box[0]],
                   [min_box[1], max_box[1]],
                   [min_box[2], max_box[2]]]

    small_boundaries = [ [min_box[0], max_box[0]],
                         [l - sph_rad, l + sph_rad],
                         [min_box[2], max_box[2]]]

    return boundaries, small_boundaries, int_spheres, sph_rad

def section(iteration_level = 3, y = 0, rotation = False, seed = None, N_points=2500):
    """
    This function returns the image resulting from the virtual section of a
    ramification.

    Parameters:
    iteration_level = The # of iterative biforcation made in the ramification
    N               = Parameter that regulate the sampling density
    """
    start = time.time()

    Pancreas = createTree(iter = iteration_level, rotation = rotation, seed = seed) #Ramification object
    boundaries, small_boundaries, int_spheres, sph_rad = _box_and_spheres(Pancreas)

    l = Pancreas[0].l

    #Defining the problem for a low discrepancy sampling inside 'boundaries'
    problem = {'num_vars': 3,
               'names': ['x', 'y', 'z'],
               'bounds': boundaries}

    N = 10
    vor_points = []

    while len(vor_points) < N_points:
        N = N*2
        vor_points = saltelli.sample(problem, N)
        vor_points = vor_points[(vor_points[:,1] > small_boundaries[1][0]) & (vor_points[:,1] < small_boundaries[1][1])]
        vor_points = vor_points[:N_points]
    vor = Voronoi(vor_points) #Creating the tassellation


    #Cropping the regions that lies outside the boundaries
    cropped_reg = [ reg for reg in vor.regions if _inside_boundaries(vor, reg, small_boundaries)]

    start_distances = time.time()
    tree = cKDTree(vor.vertices)

    #Vertices that lies inside the spheres
    inside_indexes = tree.query_ball_point(int_spheres, sph_rad)

    vertices_truth = np.zeros(len(vor.vertices)) #Array where to store identities
    for ind in inside_indexes:
        vertices_truth[ind] = 1
    vertices_truth = vertices_truth.astype(int) #0: Outside, 1: Inside

    region_id = np.zeros(len(cropped_reg))
    for n, reg in enumerate(cropped_reg):
        region_id[n] = any(vertices_truth[reg]) + all(vertices_truth[reg]) # 0:Outside, 1:Partially, 2:Inside
    region_id = region_id.astype(int)

    start_drawing = time.time()
    #DRAWING SETINGS
    fig = plt.figure(figsize=(12,8))
    colors = ['w', 'c', 'r']

    plt.axes().set_facecolor("grey")
    for n, reg in enumerate(cropped_reg):
        ind_abo = [ ver for ver in vor.vertices[reg] if ver[1] > l ]
        ind_bel = [ ver for ver in vor.vertices[reg] if ver[1] <= l ]

        couples = list(itertools.product(ind_abo,ind_bel))
        if couples:
            intersection_point = [ _plane_y_intersection(v1, v2, k = l) for v1, v2 in couples]
            intersection_point = np.asarray(intersection_point)
            drawing_hull = ConvexHull(intersection_point[:,[0,2]]).vertices
            plt.gca().fill(intersection_point[drawing_hull,0],
                           intersection_point[drawing_hull,2],
                           colors[region_id[n]])
    end = time.time()
    times = {'N'         : N_points,
             'iter_lev'  : iteration_level,
             'Voronoi'   : start_distances - start,
             'Distances' : start_drawing -  start_distances,
             'Drawing'   : end - start_drawing }

    return fig, times



#%%-----------------------------------------------------------------------------
section(iteration_level = 5, y = 2, rotation = False, seed = None, N_points=10000)

















# if __name__ == '__main__':
#     time_measures = []

    # for N_points in np.arange(25000, 2700, 100):
    #     for n in range(1):
    #         fig, times = section(rotation = True, z=0, seed = 32, N_points=N_points)
    #         tic = time.time()
    #         dpi = 100
    #         fig.savefig(f'Times\\Images\\section_N_{N_points}({n}).png', bbox_inches='tight', dpi=dpi)
    #         toc = time.time()
    #         times.update({f'Saving figure': toc-tic,
    #                       'dpi': dpi})
    #         time_measures.append(times)
    #         plt.close(fig)
    #     time_df = pd.DataFrame(time_measures)
    #     # time_df.to_csv('Times\\time_measures.csv')
    #     print(f'N = {N_points}')




# dpi = 100
# fig.savefig(f'section_N_{N_points}.png', bbox_inches='tight', dpi=dpi)
