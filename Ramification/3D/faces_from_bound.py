import itertools
from operator import itemgetter


max_box = np.max(np.asarray(centers), axis = 0) + tree_rad*2 #Padding proportional to spheres'radius
min_box = np.min(np.asarray(centers), axis = 0) - tree_rad*2
bounds = [ [min_box[i], max_box[i]] for i in range(3)]

#Finding faces of a Box
box_vertex = [el for el in itertools.product(*bounds)]
face1 = list(itemgetter(0,1,2,3)(box_vertex))
face2 = list(itemgetter(4,5,6,7)(box_vertex))
face3 = list(itemgetter(0,1,4,5)(box_vertex))
face4 = list(itemgetter(2,3,6,7)(box_vertex))
face5 = list(itemgetter(0,2,4,6)(box_vertex))
face6 = list(itemgetter(1,3,5,7)(box_vertex))
