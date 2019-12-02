from sympy import *
init_printing(use_unicode=True) #Pretty Printing

x,y,z = symbols('x y z') #Introducing main symbols

xc, yc, zc, r = (1, 2, 2, 5) #Specific variables for sphere
unit_sphere = Eq((x - xc)**2 + (y -yc)**2 + (z - zc)**2, r**2) #The actual sphere
simplify(unit_sphere.subs(z,0)) #The instersection with the plane

#%%
