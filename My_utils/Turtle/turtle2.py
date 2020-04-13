#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 13:15:20 2019

@author: alessandro
"""

import turtle

wn=turtle.Screen()
wn.bgcolor("lightgreen")
tess = turtle.Turtle() # Create tess and set some attributes
tess.color("hotpink")
tess.pensize(5)
alex = turtle.Turtle() # Create alex

tess.forward(80) # Make tess draw equilateral triangle
tess.left(120)
tess.forward(80)
tess.left(120)
tess.forward(80)
tess.left(120) # Complete the triangle

alex = turtle.forward(100)

turtle.done()
turtle.bye()
