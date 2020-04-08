import numpy as np
import itertools

from scipy.spatial import Voronoi, ConvexHull, cKDTree

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

spheres = [0,0,0]
sph_rad = 5

bounds = [[0,10],[0,10],[0,10]]

n_sample = 20
coords = np.zeros((3,n_sample))
coords[0] = np.linspace(bounds[0][0], bounds[0][1], n_sample)
coords[1] = np.linspace(bounds[1][0], bounds[1][1], n_sample)
coords[2] = np.linspace(bounds[2][0], bounds[2][1], n_sample)
reg_points = np.asarray([ pt for pt in itertools.product(coords[0], coords[1], coords[2])])
vor = Voronoi(reg_points)
len(vor.regions)

#Cropping Voronoi
crop_reg = [ reg for reg in vor.regions if inside_bounds(reg, bounds)]
len(crop_reg)

crop_ver = set()
for reg in crop_reg:
    crop_ver |= set(reg)
crop_ver = list(crop_ver)

#-------------------------------------------------------------------------------
#%% IDENTITY INSPECTION
tree = cKDTree(vor.vertices)
tree.query_ball_point(spheres, sph_rad)

 #All the spheres have the same radius
inside_indexes = tree.query_ball_point(spheres, sph_rad)
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

np.argwhere(region_id == 2)
