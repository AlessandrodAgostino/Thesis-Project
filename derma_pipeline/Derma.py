from noise import pnoise2
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

x_max = 1
y_max = 1
n_pts = 200
scale = 7 #Delicate Parameter
offset = 10 #Acts like a seed
ED_height = 100
EK_height = 20
E_thick = 100
vec_noise = np.vectorize(pnoise2)

x_val  = np.linspace(offset, x_max + offset, n_pts)
y_val  = np.linspace(offset, y_max + offset, n_pts)
xx, yy = np.meshgrid(x_val, y_val, sparse = False)
zz = vec_noise(xx, yy) * EK_height + E_thick

xx = (xx - offset) * scale + offset
yy = (yy - offset) * scale + offset
zz1 = vec_noise(xx, yy) * ED_height

EK_bound = (xx, yy, zz)
ED_bound = (xx, yy, zz1)

K_bound = None #TODO: Future iplementation

max_z = 2 * EK_height + E_thick
min_z = -ED_height*2

boundaries = np.array([[offset, offset + x_max],
                       [offset, offset + y_max],
                       [min_z , max_z       ]])

# return boundaries, ED_bound, EK_bound, K_bound

# fig = plt.figure(figsize=(15,10))
# ax = fig.gca(projection='3d')
# surf = ax.plot_trisurf(*[cc.flatten() for cc in ED_bound], label='Epid - Derm')
# surf1 = ax.plot_trisurf(*[cc.flatten() for cc in EK_bound], label='Kera - Epid')
# surf._facecolors2d=surf._facecolors3d
# surf._edgecolors2d=surf._edgecolors3d
# surf1._facecolors2d=surf1._facecolors3d
# surf1._edgecolors2d=surf1._edgecolors3d
# leg = ax.legend()
# ax.set_title("Perlin Noise", fontsize=20)
np.unique(xx)[1] - np.unique(xx)[0]
np.unique(yy)[1] - np.unique(yy)[0]
