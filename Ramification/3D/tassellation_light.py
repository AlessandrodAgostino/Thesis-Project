import numpy as np

from scipy.spatial import Voronoi, ConvexHull
from SALib.sample import saltelli

from D3branch import *

#RUN THIS SCRIPT FORM TERMINAL - it'll open a browser with Vpython GUI

#------------------------------------------------------------------------------#
#%% VORONOI PREPARATION
#Create the ramification
tree = createTree(iter = 3)

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
max_box = np.max(np.asarray(centers), axis = 0) + tree_rad*3 #Padding proportional to spheres'radius
min_box = np.min(np.asarray(centers), axis = 0) - tree_rad*3
bounds = [ [min_box[i], max_box[i]] for i in range(3)]

#SAMPLED POINTS FOR VORONOI
#Low discrepancy sampling of the box:
problem = {'num_vars': 3,
           'names': ['x', 'y', 'z'],
           'bounds': bounds}

N = 100 #No Idea how this parameter works precisely
vor_points = saltelli.sample(problem, N)

#------------------------------------------------------------------------------#
#%% VORONOI CREATION
vor = Voronoi(vor_points)

def inside_bounds(reg, bounds):
    """
    This function check if a region, hence a list of vertices, lies inside some 'boundaries'.
    """
    if -1 in reg: return False
    if len(reg)>1:
        reg = np.asarray([vor.vertices[pt] for pt in reg])
        if all(reg[:,0] > bounds[0][0]) and all(reg[:,0] < bounds[0][1]):
            if all(reg[:,1] > bounds[1][0]) and all(reg[:,1] < bounds[1][1]):
                if all(reg[:,2] > bounds[2][0]) and all(reg[:,2] < bounds[2][1]):
                    return True
        return False
    else: return False #RANDOM

crop_reg = [ reg for reg in vor.regions if inside_bounds(reg, bounds)]

#------------------------------------------------------------------------------#
#%% DRAWING METHODS:
#Creating the right canvas
scene = canvas(width=1500, height=900, center=vector(5,5,0))

#Drawing a Voronoi Tassels and their volumes if they're finite
for n,reg in enumerate(crop_reg):
    region = [vector(*vor.vertices[ver]) for ver in reg]
    drawPoints(region, color = color.hsv_to_rgb(vector (0.5,1,0.8)))

    conv_hull= ConvexHull([vor.vertices[ver] for ver in reg])
    simpl = []
    for sim in conv_hull.simplices:
        pts = [ conv_hull.points[pt] for pt in sim]
        simpl.append( triangle( vs=[vertex( pos = vector(*ver),
                                            color = color.hsv_to_rgb(vector(0.5,1,0.8)),
                                            opacity = 0.2) for ver in pts]))
print(f'The regions are:{len(crop_reg)}')

drawListBranch(tree)
drawSphereFreeEnds(tree)
