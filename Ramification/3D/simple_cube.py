import numpy as np
import itertools

from scipy.spatial import Voronoi, ConvexHull, cKDTree
from vpython import cylinder, sphere, vector, color, canvas, triangle, vertex, arrow, label


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

vor = Voronoi(reg_points) #Creating the tassellation
crop_reg = [ reg for reg in vor.regions if inside_bounds(vor, reg, bounds)]
sel_reg = crop_reg[0]
#%%-----------------------------------------------------------------------------
#The plane will be z = 7
#I'll do everything by hand then automatize it in a second moment

def plane_z_7(x,y,z):
    return (z - 7) > 0

pt_id = [ plane_z_7(*ver) for ver in vor.vertices[sel_reg] ]
ind_abo = [ n for n,ver in enumerate(vor.vertices[sel_reg]) if plane_z_7(*ver)]
ind_bel = [ n for n,ver in enumerate(vor.vertices[sel_reg]) if not plane_z_7(*ver)]

#finding every line between two points of different class
couples = list(itertools.product(ind_abo,ind_bel))

#%%-----------------------------------------------------------------------------
scene     = canvas(width=1500, height=900, center=vector(5,5,0))
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

#drawing lines between couples
for couple in couples:
    v1 = vector(*vor.vertices[sel_reg][couple[0]])
    v2 = vector(*vor.vertices[sel_reg][couple[1]])
    Figures.append(cylinder(pos = v1,
                            axis = v2 - v1,
                            radius = 0.01))


for pt in container_points:
    Figures.append(sphere(pos = vector(*pt),
                          radius = max_coord/50,
                          color = red,
                          opacity = 0.5))
