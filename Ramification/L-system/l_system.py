#@author: Alessandro d'Agostino
import cmath
import numpy as np
from numpy import random
import matplotlib.pyplot as plt
#from shapely.geometry import LineString

from branch import Branch

class L_System:
    """
    This class is intended as a set of 'Branch' objects that compose a single L-system drawing.

    - branches: the list containing all the branches
    - n_iter: the number of iterations performed on the current drawing
    """
    def __init__(self, starting_branch = None):
        self.n_iter = 0
        if starting_branch is not None:
            self.branches = [starting_branch]
        else:
            self.branches = []

    def start(self, branch = Branch()):
        """
        Method that initialize/overwrite the previous drawing.
        The default starting branch is the tail = (0,0), head = (0,1), iter_lev = 0, origin = None]
        """
        self.n_iter = 0
        self.branches.append(branch)

    def iteration(self, rule = [( + 85/180*np.pi, 1.5), (-85/180*np.pi, 1.5)], noise = True):
        """
        This method is what differentiate an L-system form another. The necessary rules for the costruction are:
            -1) number of child branch (n = 2,3)
            -2) angle deviation of each child branch (+ alpha, - alpha)
            -3) increasing/decreasing length ratio (l1 = l0 * R).
            -4) Presence of noise on the angle in the ramification. The noise it's not already adjustable.
                It's greater and greater going on with the iterations (proportional to iter_lev).

        The current branch is used to create the new branch(es). The first three rules could be all condesed
        in a list like this one:    [(+alpha, R), (-alpha, R)].
        This rule create 2 new branches (list lenght) with the reported angular deviation al lenght ratio.
        The origin and generation relationship shuold be update after every application of the rules.
        """
        np.random.seed = 30
        ang_noise = 10/180*np.pi #Fixed anguar noise to 5°
        if self.branches:
            for br in self.branches:
                if br.iter_lev == self.n_iter:
                    l0, theta0 =cmath.polar(br.head - br.tail)
                    for gen in rule:
                        l1 = l0/gen[1]
                        theta1 = theta0 + gen[0]
                        theta1 = theta1 + noise*random.uniform(-br.iter_lev*ang_noise, br.iter_lev*ang_noise)  #Adding noise
                        travel = cmath.rect(l1, theta1)
                        new_branch = Branch(tail = br.head,
                                            head = br.head + travel,
                                            iter_lev = self.n_iter + 1,
                                            origin = br)
                        br.generate.append(new_branch)
                        self.branches.append(new_branch)
                        #print("new branch: {0:.2f} to {0:.2f}".format(new_branch.tail, new_branch.head))

            self.n_iter += 1 #updating the n. of iteration
        else:
            raise Exception('Start your L_system before iteration.')

    def multiple_iterations(self, n, rule = [( + 85/180*np.pi, 1.5), (-85/180*np.pi, 1.5)], noise = True):
        """
        This method simply execute n times the itaration
        """
        if n<=10:
            for _ in range(n):
                self.iteration(rule)
        else: raise Exception('n is too high. It may take too much time')

    def _lies_between(self, P1, P2, P):
        """
        This function check if the point P lies in the plane-strip defined by the segment.
        The strip is delimited by the two lines perpendicular to the line crossing P1 and P2,
        that cross their turn P1 and P2.
        """
        #different y
        if (P2[1] - P1[1]):
            #different x
            if (P2[0] - P1[0]):
                #Probabily it's this step that doesn't work
                a = (P2[0] - P1[0]) / (P2[1] - P1[1])
                b = 1
                c1 = -P1[1] + P1[0] * (P1[0] - P2[0]) / (P2[1] - P1[1])
                c2 = -P2[1] + P2[0] * (P1[0] - P2[0]) / (P2[1] - P1[1])

                if min(-c1, -c2) < (P[0]*a + P[1]*b) < max(-c1, -c2):
                    return True
                else: return False
            #diff y, same x
            elif min(P1[1], P2[1]) < P[1] < max(P1[1], P2[1]):
                return True
            else: return False
        #same y, diff x
        elif (P2[0] - P1[0]):
            if min(P1[0], P2[0]) < P[0] < max(P1[0], P2[0]):
                return True
            else: return False
        #same y, same x
        else: raise Exception("The two points should be different.")

    def _line_point_dist(self, P1, P2, P):
        """
        If the point P lies between the plane-strip identified by the P1 P2 segment then the
        function return the distance between the line identified from P1 and P2 and the point P.
        Otherwise the function returns the minimum between the distances from P1 and P2.
        """
        if self._lies_between(P1, P2, P):
            abs = np.abs(((P2[1]-P1[1])*P[0] - (P2[0]-P1[0])*P[1] + P2[0]*P1[1] - P2[1]*P1[0]))
            den = np.sqrt(np.square(P2[1]-P1[1]) + np.square(P2[0]-P1[0]))
            return abs/den
        else:
            dist1 = np.sqrt(np.square(P[1]-P1[1]) + np.square(P[0]-P1[0]))
            dist2 = np.sqrt(np.square(P[1]-P2[1]) + np.square(P[0]-P2[0]))
            return min(dist1,dist2)

    def _spline_point_dist(self, sP, P):
        """
        This function compute the distance between a point and a splice, described by a list of points.
        """
        if len(sP)>2:
            dist = np.min([self._line_point_dist(P1,P2,P) for P1, P2 in zip(sP[:,:],sP[1:,:])])
        else:
            dist = self._line_point_dist(*sP, P)
        return dist

    def _max_radius(self,br):
        """
        Looks for the maximum radius available for the circle on a certain free end branch.
        This method is called by the 'draw' method when the 'circle' parameter is equal to 'smart'.
        The method check the distance between the free end and all the segment in
        the DIRECT dinasty of that branch.
        """
        dinasty = [] #Container for the branch's dinasty
        parent = br.origin
        while parent is not None:
            dinasty.append(parent)
            parent = parent.origin

        dinasty_coord = np.array([(br.tail.real,br.tail.imag)] + [(bra.tail.real,bra.tail.imag) for bra in dinasty])
        free_end = np.array([br.head.real, br.head.imag])
        return self._spline_point_dist(dinasty_coord, free_end) #Max distance from the free end


    def draw(self, circle = 'smart', **kwargs):
        """
        The method draw the entire L-System.
        It's possible to draw also circles at every free end of the system in
        different ways depending on the parameter 'circle':
            - 'no_circles': Any circle is drawn.
            - 'std': circles with radius equale to le free branch module are drawn in every free branch.
            - 'smart': circle with variable radius are drawn in order to avoid overlappings.

        Free end branches as 'Branch.generate == []'.
        """
        fig = plt.figure(figsize=(12,10))
        for br in self.branches:
            br.draw(**kwargs)
            if circle == 'std':
                if br.generate == []: #Free end
                    radius = abs(br.head-br.tail)
                    plt.gca().add_artist(plt.Circle((br.head.real,br.head.imag),radius ))
            elif circle == 'smart':
                if br.generate == []: #Free end
                    radius = self._max_radius(br)
                    permitted_radius = abs(br.head-br.tail)
                    if radius > permitted_radius:
                        plt.gca().add_artist(plt.Circle((br.head.real,br.head.imag),radius , color=kwargs['c']))
                    else: plt.gca().add_artist(plt.Circle((br.head.real,br.head.imag),radius , color='b'))

            elif circle == 'no_circles':
                pass #no circle is drawn
            else: raise Exception("Invalid value for 'circle' parameter")
        plt.show()


#%%
br10= Branch()
ls = L_System()
ls.start(br10)
ls.multiple_iterations(7)
# ls.draw(c = 'b', circle = 'std')
ls.draw(c = 'r', circle = 'smart')
#%%
br = ls.branches[150]
dinasty = [] #Container for the branch's dinasty
parent = br.origin
while parent is not None:
    dinasty.append(parent)
    parent = parent.origin

dinasty_coord = np.array([(br.tail.real,br.tail.imag)] + [(bra.tail.real,bra.tail.imag) for bra in dinasty])
free_end = np.array([br.head.real, br.head.imag])

sP = dinasty_coord
P = free_end

#%%
for i in range(len(sP)-1):
    print("distance from the segment = {}".format(ls._line_point_dist(sP[i],sP[i+1],P)))
    print("does it lie? {}".format(ls._lies_between(sP[i],sP[i+1],P)))
    plt.plot(sP[:,0], sP[:,1])
    plt.scatter(*free_end, c='red')
    plt.scatter(*sP[i], c='green')
    plt.scatter(*sP[i+1], c='green')
