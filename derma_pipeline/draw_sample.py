import numpy as np
import itertools

from scipy.spatial import Voronoi, ConvexHull, cKDTree, Delaunay
from SALib.sample import saltelli
from vpython import arrow, label

from ../pipeline import section
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

def draw_axis(Figures, max_coord):
    Figures.append( arrow( pos=vector(0,0,0), axis=vector(0,0,max_coord), shaftwidth=0.1))
    Figures.append( label( pos=vector(0,0,max_coord/2), text='Z' ))   #Z axis

    Figures.append( arrow( pos=vector(0,0,0), axis=vector(0,max_coord,0), shaftwidth=0.1))
    Figures.append( label( pos=vector(0,max_coord/2,0), text='Y' ))   #Y axis

    Figures.append( arrow( pos=vector(0,0,0), axis=vector(max_coord,0,0), shaftwidth=0.1))
    Figures.append( label( pos=vector(max_coord/2,0,0), text='X' ))   #X axis

def plane_z_intersection(p1, p2, z=0):
    x = p1[0] + (0 - p1[2]) / (p2[2] - p1[2]) * (p2[0] - p1[0])
    y = p1[1] + (0 - p1[2]) / (p2[2] - p1[2]) * (p2[1] - p1[1])
    return np.asarray((x,y,0))

#-------------------------------------------------------------------------------
#%% VORONOI - PREPARATION
Pancreas = createTree(iter = 2, rotation = False, seed = 49) #Ramification object

vor_points.shape


#-------------------------------------------------------------------------------
#%% VORONOI - CREATION
vor = Voronoi(vor_points) #Creating the tassellation

#Cropping the regions that lies outside the boundaries
crop_reg = [ reg for reg in vor.regions if inside_bounds(vor, reg, bounds)]
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

region_id = np.zeros(len(vor.regions))

#-------------------------------------------------------------------------------
#%% DRAWING METHODS:
"""
TIPICALLY A FILTER ON WHAT REGION TO DRAW IS NEEDED, OTHERWISE OR THE DRAWING
BECOMES MESSY OR VPYTHON DOESN'T LOAD THAT MANY OBJECTS.
The filtering is done by deciding the color of the region to draw.
Should be implemented in a more sensitive way.
"""
#Preliminary graphic settings
scene     = canvas(width=800, height=600, center=vector(5,5,0), background=color.gray(0.6))
turquoise = color.hsv_to_rgb(vector(0.5,1,0.8))
red       = color.red #Some colors
white     = color.white
black     = color.black
colors    = [black, white, turquoise,  red, black] #[Outside, Partially, Inside, Out of Boundaries]
Figures   =  [] #List to which append all the drawings

draw_axis(Figures, 10)

#Drawing Vor points depending on the color
for n,ver in enumerate(vor.points):
    if colors[vertices_truth[n]] in []: #[white, turquoise] or [orange] for nothing
        Figures.append(sphere(pos     = vector(*ver),
                              radius  = sph_rad/50,
                              opacity = 0.4,
                              color   = colors[nuclei_id[n]])) #warning

#Drawing a Voronoi Tassels and their volumes if they're finite
for n,reg in enumerate(vor.regions):
    if colors[region_id[n]] in [turquoise]: #[red, turquoise] or [orange] for nothing
        conv_hull= ConvexHull([vor.vertices[ver] for ver in reg])
        simpl = []
        for sim in conv_hull.simplices:
            pts = [conv_hull.points[pt] for pt in sim]
            simpl.append( triangle( vs=[vertex( pos     = vector(*ver),
                                                color   = colors[region_id[n]],
                                                opacity = 0.2) for ver in pts]))
