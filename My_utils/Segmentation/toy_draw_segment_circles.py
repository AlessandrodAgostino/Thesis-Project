import numpy as np
import time
from scipy.ndimage import label
import matplotlib.pyplot as plt
from shapely.ops import polygonize,unary_union
from shapely.geometry import LineString, MultiPolygon, MultiPoint, Point, Polygon
from scipy.spatial import Voronoi
from SALib.sample import saltelli
from skimage import io, morphology, img_as_uint

#structure elements for recognizing features
s = [[1,1,1],
     [1,1,1],
     [1,1,1]]

def select(pix, bot = 0, top = 255):
    if pix >= bot and pix <= top: return 1
    else: return 0
v_select = np.vectorize(select)

###DRAWING CIRCLES
circles_id = {'c1': {'r': 0.5, 'c': (2,3)},
              'b1': {'r': 1, 'c': (2,3)},
              'c2': {'r': 0.5, 'c': (7,6)},
              'b2': {'r': 1, 'c': (7,6)}}

circles = [] #list of Polygon containing circles
n_points = 30 #Circumference's Resolution
thetas = np.linspace(0,2*np.pi,n_points)
for num, c_id in circles_id.items():
    x = np.cos(thetas)*c_id['r'] + c_id['c'][0]
    y = np.sin(thetas)*c_id['r'] + c_id['c'][1]
    circle_point = np.stack((x,y)).T
    circles.append(Polygon(circle_point))

#Plotting
fig = plt.figure()
for circ in circles:
    plt.fill(*circ.exterior.xy, 'k', alpha = 0.7)
plt.ylim(0,10)
plt.xlim(0,10)
plt.axis('off')
fig.savefig('small_black_circles.png', bbox_inches='tight',dpi=50)

#%%
image = io.imread('small_black_circles.png')
image = image[:,:,0:2]
#the 0 slice is used as image
#the 1 slice is used for labeling
image[:,:,1] = np.zeros(image[:,:,1].shape)

#Selecting lumes
sel_cores = v_select(image[:,:,0], top = 30)
cores, _ = label(sel_cores, structure = s)
image[:,:,1] = image[:,:,1] + cores
np.unique(cores)
np.unique(image[:,:,1])
io.imsave('cores.png',img_as_uint(cores))

#Selecting epitelia
sel_cont = v_select(image[:,:,0], bot = 31, top = 200)
cont, _ = label(sel_cont, structure = s)
image[:,:,1] = image[:,:,1] - cont
io.imsave('count.png',img_as_uint(cont))

#Selecting background
sel_back = v_select(image[:,:,0], bot = 201)
backg, _ = label(sel_back, structure = s)
io.imsave('backg.png',img_as_uint(backg))
