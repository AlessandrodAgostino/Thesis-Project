#%%
from pyquaternion import Quaternion
from math import radians
import numpy as np
from vpython import cylinder , vector, color

v = vector(1,2,3)

l = 10
ang = 30 #angle in degrees.
x = np.array([1,0,0])
y = np.array([0,1,0])
z = np.array([0,0,1])

#(The rotations are applied from left to right)
qx_85 = Quaternion(axis=x, angle=radians(85))
qy_85 = Quaternion(axis=y, angle=radians(85))
qz_85 = Quaternion(axis=z, angle=radians(85))


st_dir = y

dir_2_1 = qz_85.rotate(st_dir)
dir_2_2 = qz_85.inverse.rotate(st_dir)

dir_3_1 = (qz_85 * qx_85).rotate(st_dir)
dir_3_2 = (qz_85 * qx_85.inverse).rotate(st_dir)
dir_3_3 = (qz_85.inverse * qx_85).rotate(st_dir)
dir_3_4 = (qz_85.inverse * qx_85.inverse).rotate(st_dir)

dir_4_1 = (qz_85 * qx_85 * qz_85).rotate(st_dir)
dir_4_2 = (qz_85 * qx_85 * qz_85.inverse).rotate(st_dir)
dir_4_3 = (qz_85 * qx_85.inverse * qz_85).rotate(st_dir)
dir_4_4 = (qz_85 * qx_85.inverse * qz_85.inverse).rotate(st_dir)


#%%

#rod = cylinder(pos=vector(0,2,1),         axis=vector(5,0,0), radius=1)
SchindlerList = []
SchindlerList.append(cylinder(pos = vector(0,0,0),
                              axis = vector(*st_dir*l),
                              radius = 1))
SchindlerList.append(cylinder(pos = vector(0,l,0),
                              axis = vector(*dir_2_1*l),
                              color = color.green,
                              radius = 1))
SchindlerList.append(cylinder(pos = vector(0,l,0),
                              axis = vector(*dir_2_2*l),
                              color = color.purple,
                              radius = 1))
#----------------------------------------------------------
SchindlerList.append(cylinder(pos = vector(0,l,0) + vector(*dir_2_1*l),
                              axis = vector(*dir_3_1*l) ,
                              color = color.yellow,
                              radius = 1))
SchindlerList.append(cylinder(pos = vector(0,l,0) + vector(*dir_2_1*l),
                              axis = vector(*dir_3_2*l),
                              color = color.magenta,
                              radius = 1))
SchindlerList.append(cylinder(pos = vector(0,l,0) + vector(*dir_2_2*l),
                              axis = vector(*dir_3_3*l),
                              color = color.orange,
                              radius = 1))
SchindlerList.append(cylinder(pos = vector(0,l,0) + vector(*dir_2_2*l),
                              axis = vector(*dir_3_4*l),
                              color = color.cyan,
                              radius = 1))
#----------------------------------------------------------
SchindlerList.append(cylinder(pos = vector(0,l,0) + vector(*dir_2_1*l) + vector(*dir_3_1*l),
                              axis = vector(*dir_4_1*l),
                              radius = 1))
SchindlerList.append(cylinder(pos = vector(0,l,0) + vector(*dir_2_1*l) + vector(*dir_3_1*l),
                              axis = vector(*dir_4_2*l),
                              radius = 1))
SchindlerList.append(cylinder(pos = vector(0,l,0) + vector(*dir_2_1*l) + vector(*dir_3_2*l),
                              axis = vector(*dir_4_3*l),
                              radius = 1))
SchindlerList.append(cylinder(pos = vector(0,l,0) + vector(*dir_2_1*l) + vector(*dir_3_2*l),
                              axis = vector(*dir_4_4*l),
                              radius = 1))
