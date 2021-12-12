# -*- coding: utf-8 -*-
"""
Created on Fri May 28 11:17:36 2021

@author: raggarwal
"""

"""
Holds the 3D rigid body dynamics expressed in quaternion space.  Supplies 
dynamics propagation functions for other modules

state: s = [x, y, z, xdot, ydot, zdot, q0, q1, q2, q3, body_omega_x, body_omega_y, body_omega_z]

where [q0, q1, q2, q3] is a unit quaternion expressing the attitude of the body
with respect to the global inertial frame

x, y, and z are in the inertial frame
omega_x, omega_y and omega_z are with respect to the body frame

"""

import Quaternion as quat;
from casadi import *

class RigidBody3D:
    
    
    
    
    def __init__(self, mass_kg,
                 Ixx_kgm2, Iyy_kgm2, Izz_kgm2,
                 thrustOffsetFromCG_m):
        
        self.m = mass_kg; # kg, mass of the vehicle
        self.r_thrust = thrustOffsetFromCG_m; # m, offset of thrust vector from CG vector.  Must be a list.
        self.g = 9.81; # m/s^2, gravitational accel
        
        self.Ixx = Ixx_kgm2;
        self.Iyy = Iyy_kgm2;
        self.Izz = Izz_kgm2;
    
    
    
    
    def ds_dt(self, s, u):
        """
        Returns the state derivative due to rigid body motion
        State s is the full 13-state vector
        Control u is the [F_t (thrust), alpha (servo 1 angle about body y), beta (servo 2 angle about body z)];
        Function does not use numpy so that Casadi can call it safely.
        
        When CASADI calls this function, the assumption is s and u are MX types
        """
        
        q       = s[6:10]; # current attitude quaternion, vector
        omega_b = s[10:13]; # current body angular velocity
        
        
        
        # Unpack controls
        F_t   = u[0]; # thurst, N
        alpha = u[1]; # servo 1 angle, rad
        beta  = u[2]; # servo 2 angle, rad
        
        
        # 1. Calculate translational position derivatives
        xdot = s[3];
        ydot = s[4];
        zdot = s[5];
        
        
        
        # 2. Calculate translational velocity derivatives
        # Compute the body force wrench
        F_wrench_body = [F_t*cos(alpha)*cos(beta),
                         F_t*cos(alpha)*sin(beta),
                        -F_t*sin(alpha)];
        
        

        # make it a quaternion form for rotation (0 real component, xyz as the imag components)
        F_wrench_body_quat_form = [0,
                                   F_t*cos(alpha)*cos(beta),
                                   F_t*cos(alpha)*sin(beta),
                                  -F_t*sin(alpha)];
        
        
        q_conj = [q[0], -q[1], -q[2], -q[3]];
        
        temp  = quat.multiply(F_wrench_body_quat_form, q); # returns a list of MX's
        F_wrench_world = quat.multiply(q_conj, temp);
        
        xdotdot = 1/self.m * F_wrench_world[1] - self.g;
        ydotdot = 1/self.m * F_wrench_world[2];
        zdotdot = 1/self.m * F_wrench_world[3];
        
        
        # 3. Calculate attitude position derivatives (quaternion)
        omega_b_quat_form = [0, omega_b[0], omega_b[1], omega_b[2]];
        
        qdot = quat.multiply(q, omega_b_quat_form);
        qdot = quat.scale(qdot, 1/2); # list of MX's
        
        
        # 4. Calculate attitude velocity derivatives (Newton-Euler)
        

        # Compute the body moment wrench with cross product (r x F).  List of MX's
        M_wrench_body = [self.r_thrust[1]*F_wrench_body[2] - self.r_thrust[2]*F_wrench_body[1],
                        -self.r_thrust[0]*F_wrench_body[2] + self.r_thrust[2]*F_wrench_body[0],
                         self.r_thrust[0]*F_wrench_body[1] - self.r_thrust[1]*F_wrench_body[0]];
    
        
        # Then, use body-frame Newton-Euler (Greenwood's Dynamics)
        #omega_b_dot = self.I_inv @ (M_wrench_body - np.cross(omega_b, self.I @ omega_b, 0, 0).reshape(3, 1));


        omega_b_dot = [1/self.Ixx*(M_wrench_body[0] - (self.Izz - self.Iyy)*omega_b[1]*omega_b[2]),
                       1/self.Iyy*(M_wrench_body[1] - (self.Ixx - self.Izz)*omega_b[0]*omega_b[2]),
                       1/self.Izz*(M_wrench_body[2] - (self.Iyy - self.Ixx)*omega_b[0]*omega_b[1])];
        
        
        
        
        
        #omega_b_dot = self.I_inv_SX @ (M_wrench_body - cross(omega_b, self.I_SX @ omega_b));

        
        # Finally, stack all of our derivatives together
        ds_dt = vertcat(xdot, ydot, zdot,
                        xdotdot, ydotdot, zdotdot,
                        qdot[0], qdot[1], qdot[2], qdot[3],
                        omega_b_dot[0], omega_b_dot[1], omega_b_dot[2]);
        
        
        return ds_dt;
        
    
    
        