#%%
from pyquaternion import Quaternion
from math import radians
import numpy as np
l = 10
ang = 30 #angle in degrees.
x = np.array([1,0,0])
y = np.array([0,1,0])
z = np.array([0,0,1])

qz_30 = Quaternion(axis=[0, 0, 1], angle=radians(30))
qz_30neg = Quaternion(axis=[0, 0, 1], angle=radians(-30))
qx_30 = Quaternion(axis=[0, 1, 0], angle=radians(30))
qx_30neg = Quaternion(axis=[0, 1, 0], angle=radians(-30))

axis = qz_30.rotate(y)
axis2 = qz_30neg.rotate(y)

axis3 = (qz_30neg * qz_30).rotate(y)
axis4 = (qz_30 * qz_30).rotate(y)

axis5 = (qz_30neg *qz_30neg ).rotate(y)
axis6 = (qz_30 * qz_30neg  ).rotate(y)

#%%
from vpython import *
#rod = cylinder(pos=vector(0,2,1),         axis=vector(5,0,0), radius=1)
SchindlerList = []
SchindlerList.append(cylinder(pos = vector(0,0,0),
                              axis = vector(0,l,0),
                              radius = 1))

SchindlerList.append(cylinder(pos = vector(0,l,0),
                              axis = vector(*l*axis),
                              radius = 0.7))
SchindlerList.append(cylinder(pos = vector(0,l,0),
                              axis = vector(*l*axis2),
                              radius = 0.7))

SchindlerList.append(cylinder(pos = vector(0,l,0) + vector(*l*axis2),
                              axis = vector(*l*axis3),
                              radius = 0.7**2))
SchindlerList.append(cylinder(pos = vector(0,l,0) + vector(*l*axis2),
                              axis = vector(*l*axis4),
                              radius = 0.7**2))
SchindlerList.append(cylinder(pos = vector(0,l,0) + vector(*l*axis),
                              axis = vector(*l*axis5),
                              radius = 0.7**2))
SchindlerList.append(cylinder(pos = vector(0,l,0) + vector(*l*axis),
                              axis = vector(*l*axis6),
                              radius = 0.7**2))
