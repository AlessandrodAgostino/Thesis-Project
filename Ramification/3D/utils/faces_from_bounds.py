import numpy as np
import itertools
from operator import itemgetter


# bounds = [ [min_box[i], max_box[i]] for i in range(3)]

#Finding faces of a Box
def faces_from_bounds(bounds):
    box_vertex = [el for el in itertools.product(*bounds)]
    face[0] = list(itemgetter(0,1,2,3)(box_vertex))
    face[1] = list(itemgetter(4,5,6,7)(box_vertex))
    face[2] = list(itemgetter(0,1,4,5)(box_vertex))
    face[3] = list(itemgetter(2,3,6,7)(box_vertex))
    face[4] = list(itemgetter(0,2,4,6)(box_vertex))
    face[5] = list(itemgetter(1,3,5,7)(box_vertex))
    return faces
