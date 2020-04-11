from noise import pnoise2
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl

x_max = 20
y_max = 20
n_pts = 200
offset = 50

x_val  = np.linspace(offset, 3*x_max + offset, 3*n_pts)
y_val  = np.linspace(0, y_max, n_pts)
xx, yy = np.meshgrid(x_val, y_val, sparse = False)

vec_noise = np.vectorize(pnoise2)
#%%
zz = vec_noise(xx, yy) * 1/6 + 1/2
zz_stack = np.dstack([zz[:,:200], zz[:,200:400], zz[:,400:]])

zz_stack[100:200, 100, :]

fig, ax = plt.subplots(figsize=(8,8))
plt.imshow(zz_stack)
#%%
zz = zz.T.reshape((200,200,3)).astype('int8')
plt.imshow(zz[:,:,0], cmap='Reds')

#%%-------------------------------------------------------------------------------
#2D Plotting
maps = ['Reds', 'Greens', 'Blues' ]
slices = [slice(0,200), slice(200 ,400), slice(400,None)]
fig, axes = plt.subplots(ncols = 3, figsize = (8*3,8))

for map, ax, sl in zip(maps, axes, slices):
    ax.contourf(xx[:,sl], yy[:,sl], zz[:,sl], cmap = plt.get_cmap(map))
    ax.get_yaxis().set_ticks([])
    ax.get_xaxis().set_ticks([])
    ax.set_title(f'{map[:-1]} Noise', fontsize=20)
fig.tight_layout()
#%%-----------------------------------------------------------------------------
#Maybe changing the Main Palette to 'RdPu' I can just create a single square
#noise with random 'offset'

x_max = 15
y_max = 15
n_pts = 200
offset = 50

cmap = plt.get_cmap('RdPu')

cmap(0.0)[:-1]
cmap(0.5)[:-1]
cmap(1.0)[:-1]

x_val  = np.linspace(offset, x_max + offset, n_pts)
y_val  = np.linspace(offset, y_max + offset, n_pts)
xx, yy = np.meshgrid(x_val, y_val, sparse = False)
vec_noise = np.vectorize(pnoise2)
zz = vec_noise(xx, yy) * 128 + 128


fig, ax = plt.subplots(1, figsize = (8,8))
ax.contourf(xx, yy, zz, cmap = cmap)
ax.get_yaxis().set_ticks([])
ax.get_xaxis().set_ticks([])
fig.tight_layout()




#%%-----------------------------------------------------------------------------
#3D Plotting
fig = plt.figure(figsize=(12,8))
ax = fig.gca(projection='3d')
ax.plot_trisurf(xx.flatten(), yy.flatten(), zz.flatten(), label='parametric curve')
ax.set_title("Perlin Noise", fontsize=20)
