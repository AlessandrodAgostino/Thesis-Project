import itertools
import random

from sympy import Plane, Point, symbols
from sympy.utilities.lambdify import lambdify


from SALib.sample import saltelli

random.seed(30)
collinear = True
section_plane = None
x = symbols('x')
y = symbols('y')
z = symbols('z')

max_coord = 10
bounds = [[0,max_coord],[0,max_coord],[0,max_coord]]

problem = {'num_vars': 3,
           'names': ['x', 'y', 'z'],
           'bounds': bounds}

while collinear:
    try:
        #Selecting three points in the volume
        pl_pt = random.sample(list(saltelli.sample(problem, 400)), 3)
        #Plane crossing through them
        section_plane = Plane(*[Point(pt) for pt in pl_pt])
        collinear = False
    except ValueError:
        print("Collinear point chosen.")

plane_function = lambdify((x,y,z), section_plane.equation(), 'numpy')
