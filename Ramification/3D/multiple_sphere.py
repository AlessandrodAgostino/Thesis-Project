import numpy as np
import itertools

from scipy.spatial import Voronoi, ConvexHull, cKDTree
from SALib.sample import saltelli

from D3branch import *
"""
Here I want to make a preliminary case of identity attribution to voronoi cells
in a volume that contains only a single sphere.
"""

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
    else: return False

#-------------------------------------------------------------------------------
#Spheres list: (x_c, y_c, z_c, r)
spheres = [(0,0,0), (0,12,0)]
radius = 5
cords = np.linspace(-6, 6, 10)
cords_y = np.linspace(-6, 18, 20)
bounds = [[-6,6], [-6,18], [-6,6]]
points = np.asarray([ pt for pt in itertools.product(cords, cords_y, cords)])

vor = Voronoi(points)
crop_reg = [ reg for reg in vor.regions if inside_bounds(reg, bounds)]

#-------------------------------------------------------------------------------
#%%
tree = cKDTree(vor.vertices)
inside_indexes = tree.query_ball_point(spheres, 5)
vertices_truth = np.zeros((vor.vertices.shape[0]))
for ind in inside_indexes:
    vertices_truth[ind] = 1

region_id = np.zeros((len(crop_reg)))
for n,reg in enumerate(crop_reg):
    # 0: outside, 1:partially, 2:inside
    region_id[n] = any(vertices_truth[reg]) + all(vertices_truth[reg])
region_id = region_id.astype(int)

#-------------------------------------------------------------------------------
#%%DRAWING
scene = canvas(width=1500, height=900, center=vector(5,5,0))
turquoise = color.hsv_to_rgb(vector(0.5,1,0.8))
red = color.red
white = color.white
colors = [turquoise, red, white]
Figures =  []

#Central sphere
for sph in spheres:
    Figures.append(sphere(pos= vector(*sph) , radius= radius, opacity = 0.8) )

#Drawing all the vertices - TO BE MODIFIED
vec_vertices = [vector(*ver) for ver in vor.vertices ]
for n,ver in enumerate(vec_vertices):
    Figures.append( sphere( pos= ver, color = white, radius = 0.1 ))

#Drawing Voronoi Cropped Region
for n,reg in enumerate(crop_reg):
    # region = [vector(*vor.vertices[ver]) for ver in reg]
    # drawPoints(region, color = colors[region_id[n]]) #Overdrawing of points

    conv_hull= ConvexHull([vor.vertices[ver] for ver in reg])
    for sim in conv_hull.simplices:
        pts = [ conv_hull.points[pt] for pt in sim]
        Figures.append( triangle( vs=[vertex( pos = vector(*ver),
                                              color = colors[region_id[n]],
                                              opacity = 0.05) for ver in pts]))

#%%
