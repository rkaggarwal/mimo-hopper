# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 22:51:08 2021

@author: raggarwal
"""


"""
This class animates a 3D rigid body given its position, attitude quaternion
and thrust vector as a function of time.  Videos are saved to the working 
directory.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection, Line3D
from matplotlib import animation
import Quaternion as quat

class Animate3DRigidBody:
    
    
    def __init__(self, x_length, y_length, z_length):
        
        self.x_length = x_length;
        self.y_length = y_length;
        self.z_length = z_length;
        
        
        self.R = np.array([  [0, 1, 0],
                [0, 0, 1],
                [1, 0, 0]]); # rotates from a world frame where X is pointing up to a 
        # world frame where Z points up.  t' = Rt.  This is so that when plotting
        # in  matplotlib +X points up.
        
        
        
    def animate3DState(self, t_,
                       xyz_tape,
                       quaternion_tape,
                       u_tape,
                       saveAsFilename = None):
       
     
        dt = t_[1] - t_[0];
        
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=int(1/dt), bitrate = -1)
        
        fig1 = plt.figure(figsize = (6, 8))
        ax1 = fig1.add_subplot(111, projection='3d')
        ax1.view_init(elev = 10, azim = 45)

        ax1.set_xlabel("World Y [m]");
        ax1.set_ylabel("World Z [m]");
        ax1.set_zlabel("World X [m]");
        ax1.set_title("Horizon Simulation")
        
        ax1.set_xlim(-5, 5);
        ax1.set_ylim(-5, 5);
        ax1.set_zlim(5, 30);
        ax1.set_box_aspect((10, 10, 25));
        
                
        # make the panes transparent
        ax1.xaxis.set_pane_color((0, 0, 0, 0))
        ax1.yaxis.set_pane_color((0, 0, 0, 0))
        ax1.zaxis.set_pane_color((0, 0, 0, 0))
        
        alpha = 0
        ax1.xaxis._axinfo["grid"]['color'] =  (0,0,0,alpha)
        ax1.yaxis._axinfo["grid"]['color'] =  (0,0,0,alpha)
        ax1.zaxis._axinfo["grid"]['color'] =  (0,0,0,alpha)
        
        ax1.grid(False)
                
         
        def rotateVectorByQuaternion(p, q):
            """
            Rotates a vector p by quaternion q which maps the world csys to the body csys
            p_rot = q'*p*q.  NOTE: this is opposite the standard form because we're interested in mapping body -> world.
            
            Input and output datatypes are row vector np arrays
            """
            
            p_quat_form = np.insert(p, 0, 0); # make p in to a quaternion form
            p_rotated = quat.multiply(quat.conjugate(q), quat.multiply(p_quat_form, q))[1:];
            return np.array(p_rotated);
            
        
        
        
        def getThrustVectorRotationMatrix(u_curr, scaling_factor = 1):
            c1 = np.cos(0); # no roll
            c2 = np.cos(scaling_factor*u_curr[0]);
            c3 = np.cos(scaling_factor*u_curr[1]);
            
            s1 = np.sin(0); # no roll
            s2 = np.sin(scaling_factor*u_curr[0]);
            s3 = np.sin(scaling_factor*u_curr[1]);
            
            R =      np.array([[c2*c3, -c2*s3, s2],
                               [c1*s3 + c3*s1*s2, c1*c3 - s1*s2*s3, -c2*s1],
                               [s1*s3 - c1*c3*s2, c3*s1 + c1*s2*s3, c1*c2]]);
            
            return R;
        
        
        
        def getVerticesAndFaces(x_curr, y_curr, z_curr, q_curr):
            
            centroid = np.array([x_curr, y_curr, z_curr]);
            
            # vertices: first rotate vertices relative to the centroid, then translate all points by the xyz location
            vertA = self.R @ (centroid + rotateVectorByQuaternion(np.array([-self.x_length/2, -self.y_length/2, -self.z_length/2]), q_curr));
            vertB = self.R @ (centroid + rotateVectorByQuaternion(np.array([-self.x_length/2, +self.y_length/2, -self.z_length/2]), q_curr));
            vertC = self.R @ (centroid + rotateVectorByQuaternion(np.array([-self.x_length/2, +self.y_length/2, +self.z_length/2]), q_curr));
            vertD = self.R @ (centroid + rotateVectorByQuaternion(np.array([-self.x_length/2, -self.y_length/2, +self.z_length/2]), q_curr));
            
            vertE = self.R @ (centroid + rotateVectorByQuaternion(np.array([self.x_length/2, -self.y_length/2, -self.z_length/2]), q_curr));
            vertF = self.R @ (centroid + rotateVectorByQuaternion(np.array([self.x_length/2, +self.y_length/2, -self.z_length/2]), q_curr));
            vertG = self.R @ (centroid + rotateVectorByQuaternion(np.array([self.x_length/2, +self.y_length/2, +self.z_length/2]), q_curr));
            vertH = self.R @ (centroid + rotateVectorByQuaternion(np.array([self.x_length/2, -self.y_length/2, +self.z_length/2]), q_curr));

            vertices = np.array([vertA, vertB, vertC, vertD, vertE, vertF, vertG, vertH]);
            
            faces = [             [vertA, vertB, vertC, vertD],
                                  [vertE, vertF, vertG, vertH],
                                  [vertA, vertB, vertF, vertE],
                                  [vertB, vertC, vertG, vertF],
                                  [vertC, vertD, vertH, vertG],
                                  [vertA, vertD, vertH, vertE] ];
            
            return vertices, faces;


        
        def getBodyAxes(x_curr, y_curr, z_curr, q_curr):
            
            centroid = np.array([x_curr, y_curr, z_curr]);
            
            bodyX_axis = centroid + rotateVectorByQuaternion(np.array([.4, 0, 0]), q_curr);
            bodyY_axis = centroid + rotateVectorByQuaternion(np.array([0, .4, 0]), q_curr);
            bodyZ_axis = centroid + rotateVectorByQuaternion(np.array([0, 0, .4]), q_curr);
            
            return self.R @ bodyX_axis, self.R @ bodyY_axis, self.R @ bodyZ_axis;
        
        # bodyCenter_to_thrust_origin = np.array([-self.x_length/2, 0, 0]);
        # thrust_vector = np.array([-self.x_length/2-.5, 0, 0]);
        

        
        def animate(i):
            

            
            centroid = self.R @ xyz_tape[:, i];
            
            ax1.collections.clear()
            ax1.lines = []
            
            curr_vertices, curr_faces = getVerticesAndFaces(xyz_tape[0, i], xyz_tape[1, i], xyz_tape[2, i],
                                                            quaternion_tape[:, i]);
            
            bodyX_axis, bodyY_axis, bodyZ_axis = getBodyAxes(xyz_tape[0, i], xyz_tape[1, i], xyz_tape[2, i],
                                                            quaternion_tape[:, i]);
            
            faceCollection = Poly3DCollection(curr_faces, facecolors = 'k', linewidths = 1, edgecolor = 'k', alpha = .20)
            
            
            
            bodyX_axis_line = Line3D(    [centroid[0], bodyX_axis[0]],
                                         [centroid[1], bodyX_axis[1]],
                                         [centroid[2], bodyX_axis[2]], color = 'r');
            
            bodyY_axis_line = Line3D(    [centroid[0], bodyY_axis[0]],
                                         [centroid[1], bodyY_axis[1]],
                                         [centroid[2], bodyY_axis[2]], color = 'g');
            
            bodyZ_axis_line = Line3D(    [centroid[0], bodyZ_axis[0]],
                                         [centroid[1], bodyZ_axis[1]],
                                         [centroid[2], bodyZ_axis[2]], color = 'b');


            
            # thrustVector_originPoint_WorldFrame   = self.centroid_W + self.R_BtoW @ R_TV @ R_B @ bodyCenter_to_thrust_origin;
            # thrustVector_terminalPoint_WorldFrame = self.centroid_W + self.R_BtoW @ R_TV @ R_B @ thrust_vector;
            
            # thrustVector_line = Line3D(  [thrustVector_originPoint_WorldFrame[0], thrustVector_terminalPoint_WorldFrame[0]],
            #                              [thrustVector_originPoint_WorldFrame[1], thrustVector_terminalPoint_WorldFrame[1]],
            #                              [thrustVector_originPoint_WorldFrame[2], thrustVector_terminalPoint_WorldFrame[2]], linewidth = 5, color = 'orange');
            
            ax1.add_collection3d(faceCollection);
            ax1.add_line(bodyX_axis_line);
            ax1.add_line(bodyY_axis_line);
            ax1.add_line(bodyZ_axis_line);
            # ax1.add_line(thrustVector_line)


        anim = animation.FuncAnimation(fig1, func = animate,
                                       frames = len(t_),
                                       interval = 1,
                                       blit = False,
                                       repeat = True)
        
        if(saveAsFilename != None):
           anim.save(saveAsFilename, writer=writer)
           
