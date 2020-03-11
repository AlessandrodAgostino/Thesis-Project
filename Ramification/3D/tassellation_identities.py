import numpy as np
import itertools

from scipy.spatial import Voronoi, ConvexHull, cKDTree
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

#-------------------------------------------------------------------------------
#%% VORONOI - PREPARATION
Pancreas = createTree(iter = 2) #Ramification object

#Extracting free end's spheres and radius
spheres  = [] #List of spheres
sph_rad  = 0
max_iter = np.log2((len(Pancreas)+1)) - 1
for br in Pancreas:
    if br.iter_lev == max_iter:
        cent = br.pos + vector(*(br.length * br.drct))
        spheres.append(np.asarray((cent.x, cent.y, cent.z)))
        if not sph_rad: sph_rad = br.length

#Boundaries for random sampling in a volume with a padding proportional to spheres' radius
max_box = np.max(np.asarray(spheres), axis = 0) + sph_rad*2
min_box = np.min(np.asarray(spheres), axis = 0) - sph_rad*2
bounds  = [ [min_box[i], max_box[i]] for i in range(3)]

#Alternative boundaries for selective drawings
partial_bounds = [ [spheres[2][i], max_box[i]] for i in range(3)]

#Defining the problem for a low discrepancy sampling inside 'bounds'
problem = {'num_vars': 3,
           'names': ['x', 'y', 'z'],
           'bounds': bounds}

#Parameter that regulate the density of the sample
N = 100000 #SHOULD UNDERSTAND BETTER HOW EXACTLY WORKS
vor_points = saltelli.sample(problem, N) #Sampling

#-------------------------------------------------------------------------------
"""
#UTILITY
#Alternative REGULAR VORONOI in the same boundaries
n_sample = 30
coords = np.zeros((3,n_sample))
coords[0] = np.linspace(bounds[0][0], bounds[0][1], n_sample)
coords[1] = np.linspace(bounds[1][0], bounds[1][1], n_sample)
coords[2] = np.linspace(bounds[2][0], bounds[2][1], n_sample)
reg_points = np.asarray([ pt for pt in itertools.product(coords[0], coords[1], coords[2])])
"""
#-------------------------------------------------------------------------------
#%% VORONOI - CREATION
vor = Voronoi(vor_points) #Creating the tassellation

#Cropping the regions that lies outside the boundaries
crop_reg = [ reg for reg in vor.regions if inside_bounds(reg, partial_bounds)]

#Detecting all the vertices that lies in/outside the boundaries
crop_ver = set()
for reg in crop_reg:
    crop_ver |= set(reg)
out_bound_ver = set(np.arange(0,len(vor.vertices))) - crop_ver

#-------------------------------------------------------------------------------
#%% IDENTITY ASSIGNMENT
tree = cKDTree(vor.vertices) #Object for an improved distance computation

#Vertices that lies inside the spheres
inside_indexes = tree.query_ball_point(spheres, sph_rad)
vertices_truth = np.zeros((vor.vertices.shape[0])) #Array where to store identities
for ind in inside_indexes:
    vertices_truth[ind] = 1
vertices_truth[list(out_bound_ver)] = -1
vertices_truth = vertices_truth.astype(int) #0: Outside, 1: Inside, -1: Outside of Boundaries

#Regions that lies completely/partially inside the spheres
region_id = np.zeros((len(crop_reg))) #Array where to store identities
for n,reg in enumerate(crop_reg):
    region_id[n] = any(vertices_truth[reg]) + all(vertices_truth[reg])
region_id = region_id.astype(int) #0: Outside, 1:Partially, 2:Inside

#-------------------------------------------------------------------------------
#%% DRAWING METHODS:
"""
TIPICALLY A FILTER ON WHAT REGION TO DRAW IS NEEDED, OTHERWISE OR THE DRAWING
BECOMES MESSY OR VPYTHON DOESN'T LOAD THAT MANY OBJECTS.
The filtering is done by deciding the color of the region to draw.
Should be implemented in a more sensitive way.
"""
#Preliminary graphic settings
scene     = canvas(width=1500, height=900, center=vector(5,5,0))
turquoise = color.hsv_to_rgb(vector(0.5,1,0.8))
red       = color.red #Some colors
white     = color.white
orange    = color.orange
black     = color.black
colors    = [white, turquoise,  red, black] #[Outside, Partially, Inside, Out of Boundaries]
Figures   =  [] #List to which append all the drawings

drawListBranch(Pancreas) #Drawing ramification
drawSphereFreeEnds(Pancreas) #Drawing free ends' spheres

#Cropped vertices
for n,ver in enumerate(vor.vertices):
    if colors[vertices_truth[n]] == orange:
        Figures.append(sphere(pos     = vector(*ver),
                              radius  = sph_rad/50,
                              opacity = 0.4,
                              color   = colors[vertices_truth[n]]))

#Drawing a Voronoi Tassels and their volumes if they're finite
for n,reg in enumerate(crop_reg):
    if colors[region_id[n]] == turquoise:
        conv_hull= ConvexHull([vor.vertices[ver] for ver in reg])
        simpl = []
        for sim in conv_hull.simplices:
            pts = [conv_hull.points[pt] for pt in sim]
            simpl.append( triangle( vs=[vertex( pos     = vector(*ver),
                                                color   = colors[region_id[n]],
                                                opacity = 0.2) for ver in pts]))
