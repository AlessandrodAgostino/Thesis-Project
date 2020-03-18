from pyquaternion import Quaternion
from math import radians
import numpy as np
from vpython import cylinder, sphere, vector, color, canvas, triangle, vertex

#Create the Canvas into which draw

class D3Branch(object):
    """docstring for 3DBranch."""

    l = 10 #Length
    r = 1 #Radius
    a = 75 #Ramification angle
    rt = 0.7 #Ratio
    qx = Quaternion(axis=np.array([1,0,0]), angle=radians(a)) #Quaternion for rotation around x of a
    qz = Quaternion(axis=np.array([0,0,1]), angle=radians(a)) #Quaternion for rotation around z of a
    st_dir = np.array([0,1,0])

    def __init__(self,
                 pos = vector(0,0,0),
                 quat = Quaternion(axis=np.array([0,0,1]), angle=0),
                 rotation = False,
                 seed = None,
                 iter_lev = 0 ):


        if seed : np.random.seed(seed)
        if rotation:
            pt_in_box = np.random.rand(3) * 2 -1
            alpha = np.random.rand(1) * 2*np.pi
            self.quat = Quaternion( axis = pt_in_box , angle=alpha)
        else: self.quat = quat

        self.pos = pos
        self.drct = self.quat.rotate(self.st_dir)
        self.iter_lev = iter_lev
        self.length = self.l*self.rt**self.iter_lev
        self.radius = self.r*self.rt**self.iter_lev

    def biforcation_quat(self, quat):
        br1 = D3Branch(pos = self.pos + vector(*self.drct)*self.length,
                       quat = self.quat * quat,
                       iter_lev = self.iter_lev + 1)
        br2 = D3Branch(pos = self.pos + vector(*self.drct)*self.length,
                       quat = self.quat * quat.inverse,
                       iter_lev = self.iter_lev + 1)
        return [br1,br2]

def drawListBranch(tree, opacity = 0.5):
    SchindlerList = []
    for br in tree:
        SchindlerList.append(cylinder(pos = br.pos,
                                      axis = vector(*br.drct*br.length),
                                      radius = br.radius,
                                      opacity = opacity))
def drawSphereFreeEnds(tree, opacity = 0.5):
    SphereList = []
    max_iter = np.log2((len(tree)+1)) - 1
    for br in tree:
        if br.iter_lev == max_iter:
            SphereList.append(sphere(pos=br.pos + vector(*br.drct*br.length),
                                     radius=br.length,
                                     opacity = 0.5))

def drawPoints(points_list, rad = 0.1, color = vector(1,0,0)):
    PointList = []
    for pt in points_list:
        PointList.append(sphere(pos= pt , radius= rad, color = color))

def createTree(iter = 9, **kwargs):
    st_br = D3Branch(**kwargs)
    tree = [st_br]
    for it in range(iter):
        add_tree = []
        if (it%2):
            for br in tree:
                if br.iter_lev == it:
                    add_tree = add_tree + br.biforcation_quat(br.qx)

        else:
            for br in tree:
                if br.iter_lev == it:
                    add_tree = add_tree + br.biforcation_quat(br.qz)
        tree = tree + add_tree

    #Return the tree as a list of D3Branch objects
    return tree

def main():
    tree = createTree(rotation = False, seed = 30)

    #DRAWING METHODS:
    scene = canvas(width=1500, height=900, center=vector(5,5,0))
    drawListBranch(tree)
    drawSphereFreeEnds(tree)

if __name__ == '__main__':
    main()
