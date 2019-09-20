#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 11:18:43 2019

@author: alessandro
"""
from dataclasses import dataclass
import numpy as np

@dataclass
class L_Point:
    """
    Define a point with the heading direction and the distance to be covered:
        x: x coordinate
        y: y coordinate
        a: [0;360] angle that describe the heading direction of the future segment
        l: lenght of the future segment coming out of the point
    """
    x: float
    y: float
    a: float
    l: float

    def ramificate(self, D = 85, R = np.sqrt(2)):
        pass
    
#%%