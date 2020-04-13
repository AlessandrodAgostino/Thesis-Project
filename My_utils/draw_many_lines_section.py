import numpy as np
import itertools
import random

from sympy import Plane, Point

from scipy.spatial import Voronoi, ConvexHull, cKDTree, Delaunay
from vpython import cylinder, sphere, vector, color, canvas, triangle, vertex, arrow, label
from SALib.sample import saltelli

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

#%%-----------------------------------------------------------------------------
max_coord = 10
bounds = [[0,max_coord],[0,max_coord],[0,max_coord]]
n_sample = 4
coords = np.zeros((3,n_sample))
coords[0] = np.linspace(bounds[0][0], bounds[0][1], n_sample)
coords[1] = np.linspace(bounds[1][0], bounds[1][1], n_sample)
coords[2] = np.linspace(bounds[2][0], bounds[2][1], n_sample)

container_points = list(itertools.product([0,max_coord],[0,max_coord],[0,max_coord]))
reg_points = np.asarray([ pt for pt in itertools.product(coords[0], coords[1], coords[2])])
#%%-----------------------------------------------------------------------------
problem = {'num_vars': 3,
           'names': ['x', 'y', 'z'],
           'bounds': bounds}

#Parameter that regulate the density of the sample
N = 10
vor_points = saltelli.sample(problem, N) #Sampling

#Alternative Line
#vor = Voronoi(reg_points) #Creating the tassellation w/ regular grid
vor = Voronoi(vor_points) #Creating the tassellation w/sampled points

crop_reg = [ reg for reg in vor.regions if inside_bounds(vor, reg, bounds)]

#%%-----------------------------------------------------------------------------
sel_z = 7
def plane_z(x,y,z, sel_z):
    return (z - sel_z) > 0

def plane_z_intersection(p1, p2, z):
    x = p1[0] + (z - p1[2]) / (p2[2] - p1[2]) * (p2[0] - p1[0])
    y = p1[1] + (z - p1[2]) / (p2[2] - p1[2]) * (p2[1] - p1[1])
    return np.asarray((x,y,z))

#Selecting by hand a particular region
sel_reg = crop_reg[7]
pt_id = [ plane_z(*ver, sel_z) for ver in vor.vertices[sel_reg] ]
ind_abo = [ n for n,ver in enumerate(vor.vertices[sel_reg]) if plane_z(*ver, sel_z)]
ind_bel = [ n for n,ver in enumerate(vor.vertices[sel_reg]) if not plane_z(*ver, sel_z)]

#%%-----------------------------------------------------------------------------
#finding every line between two points of different class
couples = list(itertools.product(ind_abo,ind_bel))
intersection_point = []
for couple in couples:
    v1 = vor.vertices[sel_reg][couple[0]]
    v2 = vor.vertices[sel_reg][couple[1]]
    intersection_point.append(plane_z_intersection(v1, v2, sel_z))
intersection_point = np.asarray(intersection_point)
intersection_point
intersectiong_triang = Delaunay(intersection_point[:,0:2])

#%%-----------------------------------------------------------------------------
scene     = canvas(width=900, height=700, center=vector(5,5,0))
turquoise = color.hsv_to_rgb(vector(0.5,1,0.8))
red       = color.red #Some colors
white     = color.white
green     = color.green
yellow     = color.yellow

Figures   =  [] #List to which append all the drawings

draw_axis(Figures,max_coord)

#drawing points ABOVE
for ver in vor.vertices[sel_reg][ind_abo]:
    Figures.append(sphere(pos = vector(*ver),
                          radius = max_coord/50,
                          color = green,
                          opacity = 0.5))

#drawing points BELOW
for ver in vor.vertices[sel_reg][ind_bel]:
    Figures.append(sphere(pos = vector(*ver),
                          radius = max_coord/50,
                          color = yellow,
                          opacity = 0.5))

# for pt in intersection_point:
#     Figures.append(sphere(pos = vector(*pt),
#                           radius = max_coord/70,
#                           color = turquoise,
#                           opacity = 0.5))

simpl = []
for sim in intersectiong_triang.simplices:
    pts = [intersectiong_triang.points[pt] for pt in sim]
    simpl.append( triangle( vs=[vertex( pos     = vector(*ver, 7),
                                        color   = turquoise,
                                        opacity = 0.7) for ver in pts]))

#drawing lines between couples
for couple in couples:
    v1 = vector(*vor.vertices[sel_reg][couple[0]])
    v2 = vector(*vor.vertices[sel_reg][couple[1]])
    Figures.append(cylinder(pos = v1,
                            axis = v2 - v1,
                            radius = 0.01,
                            opacity  =0.5))

for pt in container_points:
    Figures.append(sphere(pos = vector(*pt),
                          radius = max_coord/50,
                          color = red,
                          opacity = 0.5))

#%%-----------------------------------------------------------------------------
