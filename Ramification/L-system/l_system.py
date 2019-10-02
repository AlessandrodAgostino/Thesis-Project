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
    """

    def __init__(self, starting_branch = None):
        self.n_iter = 0
        if starting_branch is not None:
            self.branches = [starting_branch]
        else:
            self.branches = []

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
        ang_noise = 5/180*np.pi #Fixed anguar noise to 5°
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
        for _ in range(n):
            self.iteration(rule)

    def start(self, branch = Branch()):
        """
        Method that initialize/overwrite the previous drawing.
        The default starting branch is the tail = (0,0), head = (0,1), iter_lev = 0, origin = None]
        """
        self.n_iter = 0
        self.branches.append(branch)

    def draw(self, **kwargs):
        fig = plt.figure()
        for br in self.branches:
            br.draw(**kwargs)
        plt.show()


#%%
br10= Branch()
ls = L_System()
ls.start(br10)
ls.multiple_iterations(7)
ls.draw(c = 'b')
