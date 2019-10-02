#@author: Alessandro d'Agostino
import cmath
import matplotlib.pyplot as plt

class Branch:
    """
    This class intends to implements the concept of 'branche' in a more complicated L-system drawing.
    Every branch is intended as an ordered couple of points in a 2D space(about now).
    Every branch knows from who it was genereted and who it generates.
    A branch knows if it's a ramification's end (Free extreme)
    Every branch knows its own level of iteration.
    """
    # def __init__(self,
    #              x_t = 0,
    #              y_t = 0,
    #              x_h = 0,
    #              y_h = 1,
    #              iter_lev = 0,
    #              origin = None):
    #     self.tail = complex(x_t, y_t) #Tail coordinates
    #     self.head = complex(x_h, y_h) #Head coordinates
    #     self.iter_lev = iter_lev #Iteration Level
    #     self.origin = origin #Who generates this branch, if 'None' the branch is the starting one.
    #     self.generate = [] #Who is generated by thi branch.

    def __init__(self,
                 tail = complex(0,0),
                 head = complex(0,1), #as complex numbers
                 iter_lev = 0,
                 origin = None):
        self.tail = tail #Tail coordinates
        self.head = head #Head coordinates
        self.iter_lev = iter_lev #Iteration Level
        self.origin = origin #Who generates this branch, if 'None' the branch is the starting one.
        self.generate = [] #Who is generated by thi branch.

    def draw(self, **kwargs):
        plt.plot((self.tail.real,self.head.real), (self.tail.imag, self.head.imag), **kwargs)
