import numpy as np
import itertools

from scipy.spatial import Voronoi, ConvexHull, cKDTree
from SALib.sample import saltelli

from D3branch import *

#RUN THIS SCRIPT FORM TERMINAL - it'll open a browser with Vpython GUI

def inside_bounds(reg, bounds):
    """
    This function check if a region, hence a list of vertices, lies inside some
    'boundaries'.
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

#-------------------------------------------------------------------------------
#%% VORONOI PREPARATION
#Create the ramification
Pancreas = createTree(iter = 2)

#Extracting free end's spheres
#Detecting free end's spheres radius
spheres = []
sph_rad = 0
max_iter = np.log2((len(Pancreas)+1)) - 1
for br in Pancreas:
    if br.iter_lev == max_iter:
        cent = br.pos + vector(*(br.length * br.drct))
        spheres.append(np.asarray((cent.x, cent.y, cent.z)))
        if not sph_rad: sph_rad = br.length

#Boundaries for random sampling
max_box = np.max(np.asarray(spheres), axis = 0) + sph_rad*2 #Padding proportional to spheres'radius
min_box = np.min(np.asarray(spheres), axis = 0) - sph_rad*2
bounds = [ [min_box[i], max_box[i]] for i in range(3)]
partial_bounds = [ [0, max_box[i]] for i in range(3)]

#SAMPLED POINTS FOR VORONOI
#Low discrepancy sampling of the box:
problem = {'num_vars': 3,
           'names': ['x', 'y', 'z'],
           'bounds': bounds}

N = 1000 #No Idea how this parameter works precisely
vor_points = saltelli.sample(problem, N)

#-------------------------------------------------------------------------------
#Trying with regular VORONOI
n_sample = 40
coords = np.zeros((3,n_sample))
coords[0] = np.linspace(bounds[0][0], bounds[0][1], n_sample)
coords[1] = np.linspace(bounds[1][0], bounds[1][1], n_sample)
coords[2] = np.linspace(bounds[2][0], bounds[2][1], n_sample)
reg_points = np.asarray([ pt for pt in itertools.product(coords[0], coords[1], coords[2])])

#-------------------------------------------------------------------------------
#%% VORONOI CREATION
vor = Voronoi(reg_points)

#Cropping Voronoi
crop_reg = [ reg for reg in vor.regions if inside_bounds(reg, bounds)]
crop_ver = set()
for reg in crop_reg:
    crop_ver |= set(reg)
crop_ver = list(crop_ver)

#-------------------------------------------------------------------------------
#%% IDENTITY INSPECTION
tree = cKDTree(vor.vertices[crop_ver])

 #All the spheres have the same radius
inside_indexes = tree.query_ball_point(spheres, 5)
vertices_truth = np.zeros((vor.vertices.shape[0]))
for ind in inside_indexes:
    vertices_truth[ind] = 1
vertices_truth = vertices_truth.astype(int)

region_id = np.zeros((len(crop_reg)))
for n,reg in enumerate(crop_reg):
    # 0: outside, 1:partially, 2:inside
    region_id[n] = any(vertices_truth[reg]) + all(vertices_truth[reg])
    # region_id[n] = all(vertices_truth[reg])
region_id = region_id.astype(int)

# for n,reg in enumerate(crop_reg):
#     # 0: outside, 1:partially, 2:inside
#     if all(vertices_truth[reg]): print(n)

# vertices_truth[crop_reg[2033]]

#-------------------------------------------------------------------------------
#%% DRAWING METHODS:
#Preliminary setting
scene = canvas(width=1500, height=900, center=vector(5,5,0))
turquoise = color.hsv_to_rgb(vector(0.5,1,0.8))
red = color.red
white = color.white
orange = color.orange
colors = [white, turquoise,  red]
Figures =  []

# #Cropped vertices
# for n,ver in enumerate(crop_ver):
#     if colors[vertices_truth[n]] == turquoise:
#         Figures.append(sphere(pos = vector(*vor.vertices[ver]),
#                               radius = sph_rad/50,
#                               opacity = 0.4,
#                               color = colors[vertices_truth[n]]))

#Drawing a Voronoi Tassels and their volumes if they're finite
for n,reg in enumerate(crop_reg):
    if colors[region_id[n]] == red:
        conv_hull= ConvexHull([vor.vertices[ver] for ver in reg])
        simpl = []
        for sim in conv_hull.simplices:
            pts = [ conv_hull.points[pt] for pt in sim]
            simpl.append( triangle( vs=[vertex( pos = vector(*ver),
                                                color = colors[region_id[n]],
                                                opacity = 0.2) for ver in pts]))
"""
294
1210
1347
1551
1794
2003
2270
2645
3045
"""
#DRAWING SPECIFIC REGION FOR A FURTHER ANALYSIS
# selected_reg = crop_reg[1347]
# conv_hull= ConvexHull([vor.vertices[ver] for ver in selected_reg])
# simpl = []
# for sim in conv_hull.simplices:
#     pts = [ conv_hull.points[pt] for pt in sim]
#     simpl.append( triangle( vs=[vertex( pos = vector(*ver),
#                                         color = turquoise,
#                                         opacity = 1) for ver in pts]))
#
# for pt in [vor.vertices[ver] for ver in selected_reg]:
#     Figures.append(sphere(pos = vector(*pt),
#                           radius = sph_rad/50,
#                           opacity = 0.4,
#                           color = red))
# print(vertices_truth[selected_reg])

drawListBranch(Pancreas)
drawSphereFreeEnds(Pancreas)
