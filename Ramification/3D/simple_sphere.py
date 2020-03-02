import numpy as np
import time
from sympy import Plane
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, ConvexHull
from SALib.sample import saltelli

#%% PREPARING ALL THE  VORONOI ~ 20 s
#Low discrepancy sampling of the plane:
# problem = {'num_vars': 3,
#            'names': ['x', 'y', 'z'],
#            'bounds': [[0, 10],[0, 10], [0, 10]]}
#
# low_points = saltelli.sample(problem, 3)
low_points = np.array([[2.19726562, 0.96679688, 5.18554688],
       [6.76757812, 0.96679688, 5.18554688],
       [2.19726562, 2.80273438, 5.18554688],
       [2.19726562, 0.96679688, 9.07226562],
       [2.19726562, 2.80273438, 9.07226562],
       [6.76757812, 0.96679688, 9.07226562],
       [6.76757812, 2.80273438, 5.18554688],
       [6.76757812, 2.80273438, 9.07226562],
       [7.19726562, 5.96679688, 0.18554688],
       [1.76757812, 5.96679688, 0.18554688],
       [7.19726562, 7.80273438, 0.18554688],
       [7.19726562, 5.96679688, 4.07226562],
       [7.19726562, 7.80273438, 4.07226562],
       [1.76757812, 5.96679688, 4.07226562],
       [1.76757812, 7.80273438, 0.18554688],
       [1.76757812, 7.80273438, 4.07226562],
       [9.69726562, 3.46679688, 7.68554688],
       [9.26757812, 3.46679688, 7.68554688],
       [9.69726562, 5.30273438, 7.68554688],
       [9.69726562, 3.46679688, 1.57226562],
       [9.69726562, 5.30273438, 1.57226562],
       [9.26757812, 3.46679688, 1.57226562],
       [9.26757812, 5.30273438, 7.68554688],
       [9.26757812, 5.30273438, 1.57226562]])

vor = Voronoi(low_points)

# fig = plt.figure(figsize=(10,8))
# ax = fig.add_subplot(111, projection='3d')
#
# for p in low_points:

#   ax.scatter(*p, color = 'g')
#
# for p in vor.vertices:
#   ax.scatter(*p, color = 'r')
#%%
#This region is delimited.
vor_points = [vor.points[pos] for pos in vor.regions[5]]
#this is the convex hull corresponding to the voronoi tassel
hull = ConvexHull(points=vor_points)
indexes = hull.simplices[0]
vor_points[indexes[0]]
vor_points[indexes[1]]
vor_points[indexes[2]]

final_points = [vor_points[p] for p in indexes]
#computing the plae passing true those three points
plane = Plane(*[tuple(points) for points in final_points])
plane.equation()
# tom pack
