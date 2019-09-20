import numpy as np
import time
from scipy.ndimage import label, morphology
import matplotlib.pyplot as plt
from shapely.ops import polygonize,unary_union
from shapely.geometry import LineString, MultiPolygon, MultiPoint, Point, Polygon
from scipy.spatial import Voronoi
from SALib.sample import saltelli
from skimage import io, morphology, img_as_uint, img_as_ubyte, filters, color
from skimage import util

scipy.ndimage.morphology.binary_fill_holes

#structure elements for recognizing features
s = [[1,1,1],
     [1,1,1],
     [1,1,1]]

def select(pix, bot = 0, top = 255):
    if pix >= bot and pix <= top: return 1
    else: return 0
v_select = np.vectorize(select)

#%%
####BIG IMAGE
image = io.imread('bw_nuclei_bound.png')
image.shape
temp = np.zeros((image.shape[0], image.shape[1], 5))
temp.shape
temp[:,:,0] = image[:,:,0]
image = temp
"""
0: image
1: nuclei
2: cells
3: lumes
4: boundaries
"""


#%% -------- Finding all the selections ---------
# val = list(np.linspace(0,255,num=10))
# for min, max in zip(val, val[1:]):
#     sel_image = v_select(image[:,:,0], bot=min, top=max)
#     io.imsave('selections/{:.0f}_{:.2f}_selection.png'.format(min, max),
#               img_as_ubyte((sel_image>0).astype(float)))

#%%
#Working on nuclei
nuclei = io.imread('selections/0_28.33_selection.png')
np.unique(nuclei)
nuclei = morphology.binary_erosion(nuclei)
nuclei = morphology.binary_erosion(nuclei)
nuclei = morphology.binary_dilation(nuclei)
nuclei = morphology.binary_dilation(nuclei)
nuclei = nuclei*255
np.unique(nuclei)
io.imsave('nuclei.png'.format(min, max),img_as_ubyte((nuclei>0).astype(float)))

nuclei_lab, n_nuclei = label(nuclei, structure = s)
print(n_nuclei)
image[:,:,1] = image[:,:,1] + nuclei_lab
#%%
#Finding all the cells (3% error)
cells = io.imread('selections/0_28.33_selection.png')
cells = cells - nuclei
cells = morphology.binary_dilation(cells)
cells = np.invert(cells)*255
np.unique(cells)

io.imsave('cells.png',img_as_ubyte((cells>0).astype(float)))
cells_lab, n_cells = label(cells, structure = s)
print(n_cells)
image[:,:,2] = image[:,:,2] + cells_lab

#%%
#Selecting lumes
lumes = io.imread('selections/142_170.00_selection.png')
l = 20
selem = np.resize(np.array([1]*l**2), (l, l))
lumes = filters.median(lumes, selem=selem, out=None, mask=None, shift_x=False, shift_y=False, mode='nearest', cval=0.0, behavior='ndimage')
lumes = morphology.binary_fill_holes(lumes)
io.imsave('lumes.png',img_as_ubyte((lumes>0).astype(float)))

lumes_lab, n_lumes = label(lumes, structure = s)
print(n_lumes)
image[:,:,3] = image[:,:,3] + lumes_lab

#%%
#Selecting boundaries
#TO DO : insert a control on the area of the features. soppress the too small ones
bounds = io.imread('selections/57_85.00_selection.png')
l = 20
selem = np.resize(np.array([1]*l**2), (l, l))
bounds = filters.median(bounds, selem=selem, out=None, mask=None, shift_x=False, shift_y=False, mode='nearest', cval=0.0, behavior='ndimage')
io.imsave('bounds.png',img_as_ubyte((bounds>0).astype(float)))

bounds_lab, n_bounds = label(bounds, structure = s)
print(n_bounds)
image[:,:,4] = image[:,:,4] + bounds_lab
#%%
bounds_hue = np.zeros((bounds_lab.shape))
bounds_sat = np.zeros((bounds_hue.shape))
bounds_val = np.ones((bounds_hue.shape))

hue = np.linspace(0,1, num=n_bounds+2)
hue_dict = {n : h for n,h in enumerate(hue)}

for n,h in hue_dict.items():
    bounds_hue[np.where(bounds_lab == n)]= h
    bounds_hue[np.where(lumes_lab == n)] = h

bounds_sat[np.where(bounds_lab > 0)] = 1
bounds_sat[np.where(lumes_lab > 0)] = 0.3

hsv_bounds = np.stack((bounds_hue, bounds_sat, bounds_val), axis = 2)

io.imsave('rgb_bounds.png',color.hsv2rgb(img_as_ubyte(hsv_bounds)))
#%%
#-------------------------------------------------------------------------------
