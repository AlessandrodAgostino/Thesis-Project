from pyquaternion import Quaternion
from math import radians
import numpy as np
from vpython import cylinder, sphere, vector, color, canvas, triangle, vertex, arrow, label

from D3branch import createTree, drawListBranch, drawSphereFreeEnds

from section import section

# section(iteration_level = 2, N_points = 1000)
def draw_axis(max_coord):
    Figures = []
    Figures.append( arrow( pos=vector(0,0,0), axis=vector(0,0,max_coord), shaftwidth=0.1))
    Figures.append( label( pos=vector(0,0,max_coord/2), text='Z' ))   #Z axis

    Figures.append( arrow( pos=vector(0,0,0), axis=vector(0,max_coord,0), shaftwidth=0.1))
    Figures.append( label( pos=vector(0,max_coord/2,0), text='Y' ))   #Y axis

    Figures.append( arrow( pos=vector(0,0,0), axis=vector(max_coord,0,0), shaftwidth=0.1))
    Figures.append( label( pos=vector(max_coord/2,0,0), text='X' ))   #X axis
    return Figures



#DRAWING METHODS:
scene = canvas(width=1500, height=900, center=vector(5,5,0))

draw_axis(10)

tree2 = createTree(rotation = False, delta_y = 10, seed = 30)
drawListBranch(tree2)
drawSphereFreeEnds(tree2)

sphere( pos = vector(0,0,0), color = vector(1,0,0))
