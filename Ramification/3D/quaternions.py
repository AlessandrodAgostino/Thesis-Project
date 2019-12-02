from pyquaternion import Quaternion
my_quaternion = Quaternion(axis=[1, 0, 0], angle=3.14159265)
#%%
import numpy
numpy.set_printoptions(suppress=True) # Suppress insignificant values for clarity
v = numpy.array([0., 0., 1.]) # Unit vector in the +z direction
v_prime = my_quaternion.rotate(v)
v_prime

v
#%%
q1 = Quaternion(axis=[1, 0, 0], angle=3.14159265) # Rotate 180 about X
q2 = Quaternion(axis=[0, 1, 0], angle=3.14159265 / 2) # Rotate 90 about Y
q3 = q1 * q2 # Composite rotation of q1 then q2 expressed as standard multiplication
q4 = q2 * q3
v_prime = q3.rotate(v)
v_prime
v_prime = q4.rotate(v)
v_prime
#%%
#Defining the main axis
ang = 30 #angle in degrees.
x = [1,0,0]
y = [0,1,0]
z = [0,0,1]
