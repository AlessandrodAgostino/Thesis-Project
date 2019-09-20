#Trying to implement a parametric L-system
"""
Every command is a tuple: (action, arg):
	-(F,2): Forward segment of 2 * l unities, where l is the lenght unit.
	-(+,30): Right deviation of 30 degrees.
	-([, _): Stacking an action. "]" for popping.
	-("S", _): Save position

The complete instruction is a list of tuples.
A turtle finally execute a list of instr.

"""

import turtle
import numpy as np
from numpy.random import uniform

def predicate_triang(inst_list, R):
    """Function used to evolve the instruction list"""
    ram = 3
    alpha = 125
    for inst in inst_list:
        if inst[0] == "S":
            yield ("F",1)
            yield ("[")
            yield ("+",alpha,0)
            yield ("S")
            yield ("]")
            yield ("[")
            yield ("S")
            yield ("]")
            yield ("[")
            yield ("-",alpha,0)
            yield ("S")
            yield ("]")

        elif inst[0] == "F":
            yield ("F",inst[1]*R)
            
#        elif inst[0] in ["+","-"]:
#            yield (inst[0], inst[1], inst[2]+1)
                    
        else:
            yield inst
            
def predicate(inst_list, R):
    """Function used to evolve the instruction list"""
    alpha = 85
    for inst in inst_list:
        if inst[0] == "S":
            yield ("F",1)
            yield ("[")
            yield ("+",alpha,0)
            yield ("S")
            yield ("]")
            yield ("[")
            yield ("-",alpha,0)
            yield ("S")
            yield ("]")

        elif inst[0] == "F":
            yield ("F",inst[1]*R)
            
        elif inst[0] in ["+","-"]:
            yield (inst[0], inst[1], inst[2]+1)
                    
        else:
            yield inst

def angle_noise(n_it=0):
    delta = 20 / 9
    return uniform(-delta*n_it, +delta*n_it)


def draw_l_system(tartaruga, 
			instruction_list, 
			positions,
			unit_length=5, 
			unit_angle=1,
         rep = 1):
    stack = []
    for command in instruction_list:
        tartaruga.pd()
        if command[0] in ["F"]:
            tartaruga.forward(command[1]*unit_length)
        elif command[0] == "f":
            tartaruga.pu()  # pen up - not drawing
            tartaruga.forward(command[1]*unit_length)
        elif command[0] == "+":
            tartaruga.right(command[1]*unit_angle + angle_noise(rep - command[2]))
        elif command[0] == "-":
            tartaruga.left(command[1]*unit_angle  +  angle_noise(rep - command[2]))
        elif command[0] == "S":
            positions.append(tartaruga.position())
        elif command[0] == "[":
            stack.append((tartaruga.position(), tartaruga.heading()))
        elif command[0] == "]":
            tartaruga.pu()  # pen up - not drawing
            position, heading = stack.pop()
            tartaruga.goto(position)
            tartaruga.setheading(heading)
       
#%%

u_len = 5
R = 1.5
rep = 8

instr = [("S")]
for _ in range(9):
    instr = [x for x in predicate(instr, R)]

#%%

posit = []
wn=turtle.Screen()
wn.bgcolor("lightgreen")

tess = turtle.Turtle() # Create tess and set some attributes
tess.pensize(3)
tess.color("blue")
tess.speed(0)

tess.setposition(0,-15)
tess.setheading(90)
draw_l_system(tess, instr, posit, unit_length=u_len, rep = rep)

#
#for pos in posit[:5]:
#    print(pos[0], pos[1])
#    tess.setpos(pos)
#    tess.circle(u_len)

turtle.done()
turtle.bye()  











