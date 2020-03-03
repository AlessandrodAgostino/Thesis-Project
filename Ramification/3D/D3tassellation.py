import numpy as np
import itertools

from scipy.spatial import Voronoi, ConvexHull
from SALib.sample import saltelli

from D3branch import *
#------------------------------------------------------------------------------#
#%% VORONOI PREPARATION
#Create the ramification
tree = createTree()

#Extracting free end's spheres
centers = []
for br in tree:
    cent = br.pos + vector(*(br.length * br.drct))
    centers.append(np.asarray((cent.x, cent.y, cent.z)))

#Detecting free end's spheres radius
tree_rad = 0
max_iter = np.log2((len(tree)+1)) - 1
for br in tree:
    if br.iter_lev == max_iter:
        tree_rad = br.length
        break

#Boundaries for random sampling
max_box = np.max(np.asarray(centers), axis = 0) + tree_rad*2 #Padding proportional to spheres'radius
min_box = np.min(np.asarray(centers), axis = 0) - tree_rad*2
bounds = [ [min_box[i], max_box[i]] for i in range(3)]

#SAMPLED POINTS FOR VORONOI
#Low discrepancy sampling of the box:
problem = {'num_vars': 3,
           'names': ['x', 'y', 'z'],
           'bounds': bounds}

N = 50 #No Idea how this parameter works precisely
vor_points = saltelli.sample(problem, N)
internal_points = vor_points #copy for future uses
#------------------------------------------------------------------------------#
#%% VORONOI CROPPING
#Bounds for sampling of box's faces
bounds_list = []
bounds_list.append([[min_box[0], max_box[0]],
                    [min_box[1], max_box[1]],
                    [min_box[2], min_box[2] + 0.001]])
bounds_list.append([[min_box[0], max_box[0]],
                    [min_box[1], max_box[1]],
                    [max_box[2], max_box[2] + 0.001]])

bounds_list.append([[min_box[0], max_box[0]],
                    [min_box[1], min_box[1] + 0.001],
                    [min_box[2], max_box[2]] ])
bounds_list.append([[min_box[0], max_box[0]],
                    [max_box[1], max_box[1] + 0.001],
                    [min_box[2], max_box[2]] ])

bounds_list.append([[min_box[0], min_box[0] + 0.001],
                    [min_box[1], max_box[1]],
                    [min_box[2], max_box[2]] ])
bounds_list.append([[max_box[0], max_box[0] + 0.001],
                    [min_box[1], max_box[1]],
                    [min_box[2], max_box[2]] ])

problems = [{'num_vars': 3, 'names': ['x', 'y', 'z'], 'bounds': bnd} for bnd in bounds_list]
Nf = 10
#SAMPLED POINTS FROM FACES
for prob in problems:
    face_points = saltelli.sample(prob, Nf)
    vor_points = np.concatenate((vor_points, face_points))
print(f'The points onto which build Voronoi are {vor_points.shape[0]}')
#------------------------------------------------------------------------------#
#%% CROPPED-VORONOI CREATION
#vor = Voronoi(internal_points) #UNCROPPED
vor = Voronoi(vor_points)

def inside_bounds(reg, bounds):
    """
    This function check if a region, hence a list of vertices, lies inside some 'boundaries'.
    """
    if -1 in reg: return False
    if len(reg)>1:
        reg = np.asarray([vor.vertices[pt] for pt in reg])
        if any(reg[:,0] > bounds[0][0]) and any(reg[:,0] < bounds[0][1]):
            if any(reg[:,1] > bounds[1][0]) and any(reg[:,1] < bounds[1][1]):
                if any(reg[:,2] > bounds[2][0]) and any(reg[:,2] < bounds[2][1]):
                    return True
        return False
    else: return False #RANDOM

crop_reg = [ reg for reg in vor.regions if inside_bounds(reg, bounds)]


#%% Drawing a Voronoi Tassels and their volumes if they're finite
infinite_regions = []
null_regions = []

for n,reg in enumerate(crop_reg):
    if (-1 not in reg): #Useless controls but they make a double check
        if reg:
            region = [vector(*vor.vertices[ver]) for ver in reg]
            drawPoints(region, color = color.hsv_to_rgb(vector (0.5,1,0.8)))

            conv_hull= ConvexHull([vor.vertices[ver] for ver in reg])
            simpl = []
            for sim in conv_hull.simplices:
                pts = [ conv_hull.points[pt] for pt in sim]
                #simpl.append(triangle(vs=[vertex( pos = vector(*ver), color = color.hsv_to_rgb(vector(0.5,1,0.8)), opacity = 0.4) for ver in pts]))
        else: null_regions.append(n)
    else: infinite_regions.append(n)
print(f'The regions are:{len(crop_reg)}')
print('The infinite regions are:')
print(infinite_regions)
print('The null regions are:')
print(null_regions)


#------------------------------------------------------------------------------#
#%% DRAWING METHODS:
# casting the nparray list to a vector list (Drawing necessity)

# vec_internal_points = [vector(*fp) for fp in internal_points]
# drawPoints(vec_internal_points) #Centers of Voronoi w/o sampling of faces

vec_vor_points = [vector(*vp) for vp in vor_points]
drawPoints(vec_vor_points) #Centers of Voroin

vec_vor_ver = [vector(*vvp) for vvp in vor.vertices]
# drawPoints(vec_vor_ver, color = color.hsv_to_rgb(vector (0.5,1,0.8))) #Vertex of Voronoi
drawListBranch(tree)
drawSphereFreeEnds(tree)
