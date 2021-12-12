# -*- coding: utf-8 -*-
"""
Created on Fri May 28 14:29:53 2021

@author: raggarwal
"""

"""
POINT TO POINT
"""


import numpy as np
from RigidBody3D import RigidBody3D
from HorizonPlanner import HorizonPlanner
from DataPlotter import DataPlotter
import matplotlib.pyplot as plt
from Animate3DRigidBody import Animate3DRigidBody;

plt.close('all');

mass     = 3; # kg
x_length = .5; # m
y_length = .3; # m
z_length = .3; # m

# Assume the mass is evenly distributed in the rectangular prism envelope
Ixx = 1/12*mass*(y_length**2 + z_length**2);
Iyy = 1/12*mass*(x_length**2 + z_length**2);
Izz = 1/12*mass*(x_length**2 + y_length**2);

thrust_offset_from_cg = [-x_length/2, 0, 0]; # assume thrust vector acts at bottom of vehicle


mimo = RigidBody3D(mass_kg = 3,
                   Ixx_kgm2 = Ixx, Iyy_kgm2 = Iyy, Izz_kgm2 = Izz,
                   thrustOffsetFromCG_m = thrust_offset_from_cg);


N = 100; # number of knot points
t_final = 5; # seconds
dt = t_final/N;

hp = HorizonPlanner( dynamicsModel = mimo,
                     no_of_horizon_steps = N,
                     horizon_step_length_dt = dt,
                     F_t_lims = [2*9.81, 4*9.81],
                     alpha_lims = [-np.deg2rad(30), np.deg2rad(30)],
                     beta_lims = [-np.deg2rad(30), np.deg2rad(30)]);


t_ = np.linspace(0, N*dt, N+1);

x_init = 10;

s_init = [2, 0, 0, 0, 0, 0,
          1, 0, 0, 0,
          0, 0, 0];

s_final = [5, 2, -2, 0, 0, 0,
           1, 0, 0, 0,
           0, 0, 0];


HORIZON_sol = hp.planHorizon(s_init, s_final);




# %% Plot results
 
plotter = DataPlotter(dynamicsModel = mimo,
                      horizonPlanner = hp,
                      HORIZON_sol = HORIZON_sol);

plotter.plotTelemetry();




# %% Animate Results

ani =  Animate3DRigidBody(x_length, y_length, z_length);

ani.animate3DState( t_,
                    xyz_tape = HORIZON_sol[0:3, :],
                    quaternion_tape = HORIZON_sol[6:10, :],
                    u_tape = HORIZON_sol[16:19],
                    saveAsFilename = "point_to_point.mp4");
