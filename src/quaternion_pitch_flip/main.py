# -*- coding: utf-8 -*-
"""
Created on Fri May 28 14:29:53 2021

@author: raggarwal
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
t_final = 10; # seconds
dt = t_final/N;

hp = HorizonPlanner( dynamicsModel = mimo,
                     no_of_horizon_steps = N,
                     horizon_step_length_dt = dt,
                     F_t_lims = [2*9.81, 4*9.81],
                     alpha_lims = [-np.deg2rad(30), np.deg2rad(30)],
                     beta_lims  = [-np.deg2rad(30), np.deg2rad(30)]);


t_ = np.linspace(0, N*dt, N+1);

nWaypoints = 3;
q_flip_y_axis_waypoints = np.zeros((4, nWaypoints));

q_flip_y_axis_waypoints[:, 0] = np.ravel([np.cos(0/2), np.sin(0/2)*0, np.sin(0/2)*1, np.sin(0/2)*0]);
#q_flip_y_axis_waypoints[:, 1] = np.ravel([np.cos(np.pi/2/2), np.sin(np.pi/2/2)*0, np.sin(np.pi/2/2)*1, np.sin(np.pi/2/2)*0]);
q_flip_y_axis_waypoints[:, 1] = np.ravel([np.cos(np.pi/2), np.sin(np.pi/2)*0, np.sin(np.pi/2)*1, np.sin(np.pi/2)*0]);
#q_flip_y_axis_waypoints[:, 3] = np.ravel([np.cos(3*np.pi/2/2), np.sin(3*np.pi/2/2)*0, np.sin(3*np.pi/2/2)*1, np.sin(3*np.pi/2/2)*0]);
q_flip_y_axis_waypoints[:, 2] = np.ravel([np.cos(2*np.pi/2), np.sin(2*np.pi/2)*0, np.sin(2*np.pi/2)*1, np.sin(2*np.pi/2)*0]);

waypoint_horizon_indices = np.array([0, N//2, N]); # where along the horizon are these waypoints located?


q_flip_y_axis_ref = np.zeros((4, N+1));

for i in range(N+1):
    q_flip_y_axis_ref[0, i] = np.cos(i/N*(2*np.pi)/2);
    q_flip_y_axis_ref[1, i] = np.sin(i/N*(2*np.pi)/2)*0;
    q_flip_y_axis_ref[2, i] = np.sin(i/N*(2*np.pi)/2)*1;
    q_flip_y_axis_ref[3, i] = np.sin(i/N*(2*np.pi)/2)*0;
    


fig, ax = plt.subplots(constrained_layout = True);
ax.plot(t_, q_flip_y_axis_ref[0, :], label = 'q0 (real)');
ax.plot(t_, q_flip_y_axis_ref[1, :], label = 'q1 (i)');
ax.plot(t_, q_flip_y_axis_ref[2, :], label = 'q2 (j)');
ax.plot(t_, q_flip_y_axis_ref[3, :], label = 'q3 (k)');
ax.set_xlabel("time [sec]");
ax.set_ylabel("q{}");
ax.legend(loc = 'best');
ax.set_title("Flip Maneuever Quaternion Evolution")

theta_ref = np.zeros((N+1,));
for i in range(N+1):
    theta_ref[i] = np.arcsin(2*(q_flip_y_axis_ref[0, i]*q_flip_y_axis_ref[2, i] - q_flip_y_axis_ref[3, i]*q_flip_y_axis_ref[1, i]))



x_init = 10;

s_init = [x_init, 0, 0, 0, 0, 0,
          1, 0, 0, 0,
          0, 0, 0];

s_final = [x_init, 0, 0, 0, 0, 0,
           -1, 0, 0, 0,
           0, 0, 0];


xyz_ref = np.zeros((3, N+1));
loop_rad = 5; # m

for i in range(N+1):
    xyz_ref[0, i] = (x_init + loop_rad) + loop_rad*np.sin(i/N*2*np.pi - np.pi/2);
    xyz_ref[1, i] = 0;
    xyz_ref[2, i] = loop_rad*np.cos(i/N*2*np.pi-np.pi/2);
    

HORIZON_sol = hp.planHorizon(s_init, s_final, waypoint_horizon_indices, q_flip_y_axis_waypoints, xyz_ref, q_flip_y_axis_ref);




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
                    saveAsFilename = "pitch_flip.mp4");
