#@author: Alessandro d'Agostino
import cmath
import numpy as np
from numpy import random
import matplotlib.pyplot as plt
from branch import Branch

class L_System:
    """
    This class is intended as a set of 'Branch' objects that compose a single L-system drawing.
    - branches: the list containing all the branches
    - n_iter: the number of iterations performed on the current drawing

    Usage example:
        ls = L_System()
        ls.start()
        ls.multiple_iterations(8)
        ls.draw(c = 'r', circle = 'smart')
        ls.savefig('RamificationSmartCircles.png', dpi=1000)
    """
    def __init__(self, starting_branch = None, ang_noise = 5):
        self.n_iter = 0
        self.figure, self.axes = plt.subplots(figsize=(12,10))

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
        self.branches = [branch]

    def iteration(self, rule = [( + 85/180*np.pi, 1.5), (-85/180*np.pi, 1.5)], noise = True):
        """
        This method is what differentiate an L-system form another. The necessary rules for the costruction are:
            -1) number of child branch (n = 2,3)
            -2) angle deviation of each child branch (+ alpha, - alpha)
            -3) increasing/decreasing length ratio (l1 = l0 * R).
            -4) Presence of noise on the angle in the ramification.
                [It's greater and greater going on with the iterations (proportional to iter_lev)]

        The current branch is used to create the new branch(es). The first three rules could be all condesed
        in a list, like this one:    [(+alpha, R), (-alpha, R)].
        This rule create 2 new branches (list lenght) with the reported angular deviation al lenght ratio.
        The origin and generation relationship shuold be update after every application of the rules.
        """
        if self.branches != []:
            for br in self.branches:
                if br.iter_lev == self.n_iter:
                    l0, theta0 =cmath.polar(br.head - br.tail)
                    for gen in rule:
                        ang_noise = noise*br.iter_lev*gen[0]/20
                        l1 = l0/gen[1]
                        theta1 = theta0 + gen[0]
                        theta1 = theta1 + random.uniform(-ang_noise, ang_noise)  #Adding noise
                        travel = cmath.rect(l1, theta1)
                        new_branch = Branch(tail = br.head,
                                            head = br.head + travel,
                                            iter_lev = self.n_iter + 1,
                                            origin = br)
                        br.generate.append(new_branch)
                        self.branches.append(new_branch)

            self.n_iter += 1 #updating n_iter
        else:
            raise Exception('Start your L_system before iteration.')

    def multiple_iterations(self, n, rule = [( + 85/180*np.pi, 1.5), (-85/180*np.pi, 1.5)], noise = True):
        """
        This method simply execute n times the itaration.
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
        The method manages also vertical and horizzontal segments.
        """
        #different y
        if (P2[1] - P1[1]):
            #different x
            if (P2[0] - P1[0]):
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

    def _line_point_dist(self, br, fe_br, thick = 1/6):
        """
        If the free end `fe_br` branch's head lies between the plane-strip identified by the branch `br` then the
        function returns the distance between the line identified from `br` and the point `fe_br`'s head.
        Otherwise the function returns the minimum between the distances from `br`'s extremes.
        """
        P1 = np.array([br.tail.real,br.tail.imag])
        P2 = np.array([br.head.real,br.head.imag])
        P  = np.array([fe_br.head.real,fe_br.head.imag])

        if self._lies_between(P1, P2, P):
            abs = np.abs(((P2[1]-P1[1])*P[0] - (P2[0]-P1[0])*P[1] + P2[0]*P1[1] - P2[1]*P1[0]))
            den = np.abs(br.head - br.tail)
            if thick: #If `thick` the distance should be decreased by the branch's thickness
                thickness = np.abs(br.head - br.tail) * thick
                return abs/den-thickness/2
            else:
                return abs/den
        else:
            dist1 = np.sqrt(np.square(P[1]-P1[1]) + np.square(P[0]-P1[0]))
            dist2 = np.sqrt(np.square(P[1]-P2[1]) + np.square(P[0]-P2[0]))
            return min(dist1,dist2)

    def _spline_point_dist(self, dinasty, fe_br, thick = 1/6):

        if len(dinasty) > 1:
            dist = np.min([self._line_point_dist(br, fe_br, thick) for br in dinasty])
        else:
            dist = self._line_point_dist(dinasty, fe_br, thick)
        return dist

    def _return_descent(self,br):
        """
        Generator that returns the content of `br.generate` recursively, so all the descent is yielded.
        """
        if br.generate == []:
            yield br
        else:
            for gen in br.generate:
                yield from self._return_descent(gen)

    def _find_siblings(self, br):
        """
        Returns a list with all the siblings of the fixed degree = n_iter -2
        n° siblings are all the free ends descending from the common anchestor found climbing back the tree by n steps from `br`.
        """
        sib_lev = self.n_iter - 2
        com_des = br
        if sib_lev < self.n_iter:
            for lev in range(sib_lev):
                com_des = com_des.origin
        else: raise Exception("`sib_lev` parameter's too high. Max possible value is {}".format(n_iter -1))

        siblings = [sb for sb in self._return_descent(com_des)]
        siblings.remove(br) #Removing the calling branch from siblings
        return siblings

    def _siblings_dist(self, br):
        """
        This method returns the minumum distance from any sibling.
        """
        dists = [abs(br.head - sb.head) for sb in self._find_siblings(br) ]
        return min(dists)

    def _max_radius(self,br, thick = 1/6):
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

        #Distance from the spline defined by the dinasty points.
        spl_dist = self._spline_point_dist(dinasty, br, thick)

        # TO DO: adding a list of rays, to avoid redundant computations
        # About now this method doesn't really compute in a smart way...
        #Distance from the siblings free ends.
        sib_dist = self._siblings_dist(br)
        return min(spl_dist, sib_dist/2) #The available radius is half the distance by the centers

    def _compute_size(self):
        """
        This methos returns the size
        """
        x_min = min([br.head.real for br in self.branches])
        x_max = max([br.head.real for br in self.branches])

        y_min = min([br.head.imag for br in self.branches])
        y_max = max([br.head.imag for br in self.branches])

        pad = 2*abs(self.branches[-1].head - self.branches[-1].tail)

        return ([x_min-pad, x_max+pad],[y_min-pad, y_max+pad] )

    def draw(self, circle = 'smart', thick = 1/6, **kwargs):
        """
        The method draw the entire L-System.
        It's possible to draw also circles at every free end of the system in
        different ways depending on the parameter 'circle':
            - 'no_circles': Any circle is drawn.
            - 'std': circles with radius equale to le free branch module are drawn in every free branch.
            - 'smart': circle with variable radius are drawn in order to avoid overlappings.

        Free end branches as 'Branch.generate == []'.
        """
        for br in self.branches:
            if thick:
                br.thick_draw(thick,**kwargs)
            else: br.draw(**kwargs)

            if br.generate == []: #Drawing circles at free ends
                if circle == 'std':
                    radius = abs(br.head-br.tail)
                    plt.gca().add_artist(plt.Circle((br.head.real,br.head.imag),radius ))

                elif circle == 'smart':
                    radius = self._max_radius(br, thick)
                    plt.gca().add_artist(plt.Circle((br.head.real,br.head.imag),radius , color=kwargs['c']))

                elif circle == 'no_circles':
                    pass #no circle is drawn

                else: raise Exception("Invalid value for 'circle' parameter")
        #Fixing the image size
        x_extr, y_extr = self._compute_size()
        self.axes.set_xlim(x_extr)
        self.axes.set_ylim(y_extr)
        plt.gca().axis('off')

    def savefig(self, name,**kwargs):
        self.figure.savefig(name, **kwargs)

#%%
ls = L_System()
ls.start()
alpha = 85
ls.multiple_iterations(3, rule = [( + alpha/180*np.pi, 1.5), (-alpha/180*np.pi, 1.5)], noise = True)
ls.draw(c = 'r', circle = 'smart', thick = 1/100)
