import numpy as np
import itertools
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


from scipy.spatial import Voronoi, ConvexHull, cKDTree, Delaunay
from SALib.sample import saltelli

from D3branch import *

#RUN THIS SCRIPT FORM TERMINAL - it'll open a browser with Vpython GUI

def inside_bounds(vor, reg, bounds):
    """
    This function check if a 'region' of a given 'Voronoi' tassellation lies
    inside some given 'boundaries'.

    Parameters:
    vor    = scipy.spatial.Voronoi object
    reg    = subset of vor.vertices to check
    bounds = [ [min_x, max_x], [min_y, max_y], [min_z, max_z] ]

    Return: True if 'all' the verticies lies inside the region. False otherwise
    """
    if -1 in reg: return False
    if len(reg)>1:
        reg = np.asarray([vor.vertices[pt] for pt in reg])
        if all(reg[:,0] > bounds[0][0]) and all(reg[:,0] < bounds[0][1]):
            if all(reg[:,1] > bounds[1][0]) and all(reg[:,1] < bounds[1][1]):
                if all(reg[:,2] > bounds[2][0]) and all(reg[:,2] < bounds[2][1]):
                    return True
        return False
    else: return False

def plane_z_intersection(p1, p2, z=0):
    x = p1[0] + (z - p1[2]) / (p2[2] - p1[2]) * (p2[0] - p1[0])
    y = p1[1] + (z - p1[2]) / (p2[2] - p1[2]) * (p2[1] - p1[1])
    return np.asarray((x,y,z))

#-------------------------------------------------------------------------------
#%% VORONOI - PREPARATION
#GOOD RESULTS seed = 32
Pancreas = createTree(iter = 2, rotation = True, seed = 32) #Ramification object

#Extracting free end's spheres and radius
spheres  = [] #List of spheres
sph_rad  = 0
max_iter = np.log2((len(Pancreas)+1)) - 1
for br in Pancreas:
    if br.iter_lev == max_iter:
        cent = br.pos + vector(*(br.length * br.drct))
        spheres.append(np.asarray((cent.x, cent.y, cent.z)))
        if not sph_rad: sph_rad = br.length

#Interesting spheres. Those around x,y-plane
int_spheres = [sph for sph in spheres if (sph[2] > -sph_rad) & (sph[2] < sph_rad)]

#Boundaries for random sampling in a volume with a padding proportional to spheres' radius
max_box = np.max(np.asarray(spheres), axis = 0) + sph_rad*2
min_box = np.min(np.asarray(spheres), axis = 0) - sph_rad*2
bounds  = [ [min_box[i], max_box[i]] for i in range(3)]
new_bounds = [ [min_box[0], max_box[0]],
               [min_box[1], max_box[1]],
               [-sph_rad/2
                 , sph_rad/2
                 ]]

#Defining the problem for a low discrepancy sampling inside 'bounds'
problem = {'num_vars': 3,
           'names': ['x', 'y', 'z'],
           'bounds': new_bounds}

#Parameter that regulate the sampling density
N = 7000 #SHOULD UNDERSTAND BETTER HOW EXACTLY WORKS

vor_points = saltelli.sample(problem, N) #Sampling
# vor_points = vor_points[(vor_points[:,2] < sph_rad) & (vor_points[:,2] > -sph_rad),:]
# vor_points = vor_points[(vor_points[:,2]<3) & (vor_points[:,2]> -3),:]
#-------------------------------------------------------------------------------
#%% VORONOI - CREATION
vor = Voronoi(vor_points) #Creating the tassellation

#Cropping the regions that lies outside the boundaries
crop_reg = [ reg for reg in vor.regions if inside_bounds(vor, reg, new_bounds)]

#Detecting all the vertices that lies in/outside the boundaries
crop_ver = set()
for reg in crop_reg:
    crop_ver |= set(reg)
out_bound_ver = set(np.arange(0,len(vor.vertices))) - crop_ver

#-------------------------------------------------------------------------------
#%% IDENTITY ASSIGNMENT + DRAWING
tree = cKDTree(vor.vertices) #Object for an improved distance computation


#Vertices that lies inside the spheres
inside_indexes = tree.query_ball_point(int_spheres, sph_rad)
vertices_truth = np.zeros((vor.vertices.shape[0])) #Array where to store identities
for ind in inside_indexes:
    vertices_truth[ind] = 1
vertices_truth[list(out_bound_ver)] = -1
vertices_truth = vertices_truth.astype(int) #0: Outside, 1: Inside, -1: Outside of Boundaries

region_id = np.zeros(len(vor.regions))
for n,reg in enumerate(vor.regions):
    if reg in crop_reg:
        region_id[n] = any(vertices_truth[reg]) + all(vertices_truth[reg]) +1
    else: region_id[n] = 0 # 0:Cropped, 1:Outside, 2:Partially, 3:Inside
region_id = region_id.astype(int)

#%%-----------------------------------------------------------------------------
#DRAWING SETINGS
plt.figure(figsize=(12,8))
colors = ['k', 'w', 'c', 'r', 'y', 'k']

plt.axes().set_facecolor("grey")
useless_couples = 0
for n, reg in enumerate(vor.regions):
    if region_id[n]:
        ind_abo = [ ver for ver in vor.vertices[reg] if ver[2] > 0 ]
        ind_bel = [ ver for ver in vor.vertices[reg] if ver[2] <= 0 ]

        couples = list(itertools.product(ind_abo,ind_bel))

        if couples:
            intersection_point = [ plane_z_intersection(v1, v2) for v1, v2 in couples]
            intersection_point = np.asarray(intersection_point)
            drawing_hull = ConvexHull(intersection_point[:,0:2]).vertices
            plt.gca().fill(intersection_point[drawing_hull,0],
                           intersection_point[drawing_hull,1],
                           colors[region_id[n]])

        else:
            useless_couples += 1

print(f'useless couple computations = {useless_couples}')
plt.savefig(f'specimen{N}.png', bbox_inches='tight', dpi=1000)
