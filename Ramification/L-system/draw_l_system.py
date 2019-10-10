from l_system import L_System
from branch import Branch
import numpy as np

ls = L_System()
ls.start()
ls.multiple_iterations(4, rule = [( + 40/180*np.pi, 1.5), (-40/180*np.pi, 1.5)], noise = True)
ls.draw(c = 'r', circle = 'smart')
