import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

N_points = 100
d = 3
s_0 = 0.5
x=2.0000
for i in range(10):
    x = pow(1+x,1/(d+1))
g = x

alphas = np.power(1/g, np.arange(1,d+1))
z = np.mod(np.outer(alphas,np.arange(1, N_points+1)) + s_0, 1)
z.shape
#%%
small_boundaries = np.asarray([[-1,1],[0,6],[2,4.5]])
small_boundaries
lengths = small_boundaries[:,1] - small_boundaries[:,0]
lengths

z = ((z.T * lengths) - small_boundaries[:,0]).T

#%%
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(z[0,:], z[1,:], z[2,:])
plt.show()
