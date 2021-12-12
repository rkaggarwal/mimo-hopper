# -*- coding: utf-8 -*-
"""
Created on Fri May 28 09:05:20 2021

@author: raggarwal
"""

#import numpy as np
from casadi import *


    
# %% General quaternion functions
# Inputs must be vector-types -- not quaternion class objects

def multiply(q, r):
    """
    Performs the quaternion product between q and r (q*r)
    https://www.mathworks.com/help/aeroblks/quaternionmultiplication.html
    Returns a simple list for Casadi calls.
    
    q and r are both lists of MX's
    and the return value is a list of MX's
    """
    
    qr_0 = r[0]*q[0] - r[1]*q[1] - r[2]*q[2] - r[3]*q[3];
    qr_1 = r[0]*q[1] + r[1]*q[0] - r[2]*q[3] + r[3]*q[2];
    qr_2 = r[0]*q[2] + r[1]*q[3] + r[2]*q[0] - r[3]*q[1];
    qr_3 = r[0]*q[3] - r[1]*q[2] + r[2]*q[1] + r[3]*q[0];
    
    return [qr_0, qr_1, qr_2, qr_3];
    

def scale(q, scale_factor):
    """
    Scales a quaternion element-wise by the scale factor
    """
    
    return [q[0]*scale_factor, q[1]*scale_factor, q[2]*scale_factor, q[3]*scale_factor];
    

def conjugate(q):
    """
    Returns the conjugate of this quaternion (opposite rotation)
    """
    
    return [q[0], -q[1], -q[2], -q[3]];
    