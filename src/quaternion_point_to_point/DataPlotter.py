# -*- coding: utf-8 -*-
"""
Created on Sun May 30 12:38:16 2021

@author: raggarwal
"""

import numpy as np
import matplotlib.pyplot as plt;

class DataPlotter:
    
    
    def __init__(self, dynamicsModel ,
                       horizonPlanner,
                       HORIZON_sol):
        
        self.dynamicsModel = dynamicsModel;
        self.horizonPlanner = horizonPlanner;
        self.HORIZON_sol = HORIZON_sol;
        
        self.t = np.linspace(0, horizonPlanner.horizon_length_t, horizonPlanner.N+1);
        
        self.x = HORIZON_sol[0, :];
        self.y = HORIZON_sol[1, :];
        self.z = HORIZON_sol[2, :];
        
        self.xdot = HORIZON_sol[3, :];
        self.ydot = HORIZON_sol[4, :];
        self.zdot = HORIZON_sol[5, :]; 
        
        self.q0 = HORIZON_sol[6, :];
        self.q1 = HORIZON_sol[7, :];
        self.q2 = HORIZON_sol[8, :];
        self.q3 = HORIZON_sol[9, :];
        
        self.omega_x = HORIZON_sol[10, :];
        self.omega_y = HORIZON_sol[11, :];
        self.omega_z = HORIZON_sol[12, :];
        
        self.phi = HORIZON_sol[13, :]; # roll
        self.theta = HORIZON_sol[14, :]; # pitch
        self.psi = HORIZON_sol[15, :]; # yaw
        
        self.F_t = HORIZON_sol[16, :];
        self.alpha = HORIZON_sol[17, :];
        self.beta = HORIZON_sol[18, :];
        
        
        
        
        
    def plotTelemetry(self):
        
        linewidth_thin = 1;
        linewidth_thick = 1.5;
        
        # Plot translation states & state derivs
        fig1, ax1 = plt.subplots(3, 2, figsize = (8, 5), constrained_layout = True);
        fig1.suptitle("Translation States");
        
        ax1[0, 0].plot(self.t, self.x, color = 'b', linewidth = linewidth_thin);
        ax1[1, 0].plot(self.t, self.y, color = 'b', linewidth = linewidth_thin);
        ax1[2, 0].plot(self.t, self.z, color = 'b', linewidth = linewidth_thin);
        ax1[0, 1].plot(self.t, self.xdot, color = 'b', linewidth = linewidth_thin);
        ax1[1, 1].plot(self.t, self.ydot, color = 'b', linewidth = linewidth_thin);
        ax1[2, 1].plot(self.t, self.zdot, color = 'b', linewidth = linewidth_thin);
        
        ax1[0, 0].set_ylabel("x [m]");
        ax1[1, 0].set_ylabel("y [m]");
        ax1[2, 0].set_ylabel("z [m]");
        ax1[2, 0].set_xlabel("Time [sec]");
        
        ax1[0, 1].set_ylabel("xdot [m/s]");
        ax1[1, 1].set_ylabel("ydot [m/s]");
        ax1[2, 1].set_ylabel("zdot [m/s]");
        ax1[2, 1].set_xlabel("Time [sec]");
        
        ax1[0, 0].set_ylim(bottom = 0)
        
        
        # Plot the rotation states
        fig2, ax2 = plt.subplots(3, 2, figsize = (8, 5), constrained_layout = True);
        fig2.suptitle("Rotation States");
        
        ax2[0, 0].plot(self.t, np.rad2deg(self.phi), color = 'b', linewidth = linewidth_thin);
        ax2[1, 0].plot(self.t, np.rad2deg(self.theta), color = 'b', linewidth = linewidth_thin);
        ax2[2, 0].plot(self.t, np.rad2deg(self.psi), color = 'b', linewidth = linewidth_thin);
        ax2[0, 1].plot(self.t, np.rad2deg(self.omega_x), color = 'b', linewidth = linewidth_thin);
        ax2[1, 1].plot(self.t, np.rad2deg(self.omega_y), color = 'b', linewidth = linewidth_thin);
        ax2[2, 1].plot(self.t, np.rad2deg(self.omega_z), color = 'b', linewidth = linewidth_thin);
        
        ax2[0, 0].set_ylabel("phi (roll) [deg]");
        ax2[1, 0].set_ylabel("theta (pitch) [deg]");
        ax2[2, 0].set_ylabel("psi (yaw) [deg]");
        ax2[2, 0].set_xlabel("Time [sec]");
        
        ax2[0, 1].set_ylabel("omega_x [deg/s]");
        ax2[1, 1].set_ylabel("omega_y [deg/s]");
        ax2[2, 1].set_ylabel("omega_z [deg/s]");
        ax2[2, 1].set_xlabel("Time [sec]");
        
        ax2[0, 0].set_ylim(-1.1*np.max(np.abs(np.rad2deg(self.phi))), 1.1*np.max(np.abs(np.rad2deg(self.phi))));
        ax2[1, 0].set_ylim(-1.1*np.max(np.abs(np.rad2deg(self.theta))), 1.1*np.max(np.abs(np.rad2deg(self.theta))));
        ax2[2, 0].set_ylim(-1.1*np.max(np.abs(np.rad2deg(self.psi))), 1.1*np.max(np.abs(np.rad2deg(self.psi))));
        ax2[0, 1].set_ylim(-1.1*np.max(np.abs(np.rad2deg(self.omega_x))), 1.1*np.max(np.abs(np.rad2deg(self.omega_x))));
        ax2[1, 1].set_ylim(-1.1*np.max(np.abs(np.rad2deg(self.omega_y))), 1.1*np.max(np.abs(np.rad2deg(self.omega_y))));
        ax2[2, 1].set_ylim(-1.1*np.max(np.abs(np.rad2deg(self.omega_z))), 1.1*np.max(np.abs(np.rad2deg(self.omega_z))));   
        
        
        
        # Plot the quaternion evolution
        fig4, ax4 = plt.subplots(constrained_layout = True);
        ax4.set_title("Quaternion State");
        
        ax4.plot(self.t, self.q0, label = 'q0 (real)');
        ax4.plot(self.t, self.q1, label = 'q1 (i)');
        ax4.plot(self.t, self.q2, label = 'q2 (j)');
        ax4.plot(self.t, self.q3, label = 'q3 (k)');
        ax4.legend(loc = 'best');
        ax4.set_xlabel("time [sec]");
        ax4.set_ylabel("q{}");
        
        
        
        
        
        # Plot the controls
        fig3, ax3 = plt.subplots(3, 1, figsize = (5, 5), constrained_layout = True);
        fig3.suptitle("Control Efforts");
        
        ax3[0].plot(self.t, self.F_t, color = 'r');
        ax3[1].plot(self.t, np.rad2deg(self.alpha), color = 'r');
        ax3[2].plot(self.t, np.rad2deg(self.beta), color = 'r');
        
        ax3[2].set_xlabel("Time [sec]");
        ax3[0].set_ylabel("Thrust [N]");
        ax3[1].set_ylabel("alpha [deg]");
        ax3[2].set_ylabel("beta [deg]")
        
        ax3[0].axhline(self.horizonPlanner.F_t_lims[0], color = 'k', alpha = .5, label = "F_t_min");
        ax3[0].axhline(self.horizonPlanner.F_t_lims[1], color = 'k', alpha = .5, label = "F_t_max");
        
        ax3[1].axhline(np.rad2deg(self.horizonPlanner.alpha_lims[0]), color = 'k', alpha = .5, label = "alpha_min");
        ax3[1].axhline(np.rad2deg(self.horizonPlanner.alpha_lims[1]), color = 'k', alpha = .5, label = "alpha_max");
        
        ax3[2].axhline(np.rad2deg(self.horizonPlanner.beta_lims[0]), color = 'k', alpha = .5, label = "beta_min");
        ax3[2].axhline(np.rad2deg(self.horizonPlanner.beta_lims[1]), color = 'k', alpha = .5, label = "beta_max");
        
        ax3[0].set_ylim(0, 1.25*self.horizonPlanner.F_t_lims[1]);
        ax3[1].set_ylim(1.25*np.rad2deg(self.horizonPlanner.alpha_lims[0]), 1.25*np.rad2deg(self.horizonPlanner.alpha_lims[1]));
        ax3[1].set_ylim(1.25*np.rad2deg(self.horizonPlanner.beta_lims[0]), 1.25*np.rad2deg(self.horizonPlanner.beta_lims[1]));
        
        ax3[0].legend(loc = 'best');
        ax3[1].legend(loc = 'best');
        ax3[2].legend(loc = 'best');
        
        
        
        
        # plot states and controls
        fig4, ax4 = plt.subplots(3, 3, figsize = (12, 5), constrained_layout = True);
        fig4.suptitle("Point-to-Point Trajectory")
        
        
        # First column is translation states
        ax4[0, 0].plot(self.t, self.x, color = 'b', linewidth = linewidth_thin);
        ax4[1, 0].plot(self.t, self.y, color = 'b', linewidth = linewidth_thin);
        ax4[2, 0].plot(self.t, self.z, color = 'b', linewidth = linewidth_thin);
        
        ax4[0, 0].set_ylabel("x [m]");
        ax4[1, 0].set_ylabel("y [m]");
        ax4[2, 0].set_ylabel("z [m]");
        ax4[2, 0].set_xlabel("Time [sec]");
        
        
        
        
        # Second column is rotational states
        ax4[0, 1].plot(self.t, np.rad2deg(self.phi), color = 'b', linewidth = linewidth_thin);
        ax4[1, 1].plot(self.t, np.rad2deg(self.theta), color = 'b', linewidth = linewidth_thin);
        ax4[2, 1].plot(self.t, np.rad2deg(self.psi), color = 'b', linewidth = linewidth_thin);
        
        ax4[0, 1].set_ylabel("phi (roll) [deg]");
        ax4[1, 1].set_ylabel("theta (pitch) [deg]");
        ax4[2, 1].set_ylabel("psi (yaw) [deg]");
        ax4[2, 1].set_xlabel("Time [sec]");
        



        
        # Third column is controls
        ax4[0, 2].plot(self.t, self.F_t, color = 'r');
        ax4[1, 2].plot(self.t, np.rad2deg(self.alpha), color = 'r');
        ax4[2, 2].plot(self.t, np.rad2deg(self.beta), color = 'r');
        
        ax4[2, 2].set_xlabel("Time [sec]");
        ax4[0, 2].set_ylabel("Thrust [N]");
        ax4[1, 2].set_ylabel("alpha [deg]");
        ax4[2, 2].set_ylabel("beta [deg]")
        
        ax4[0, 2].axhline(self.horizonPlanner.F_t_lims[0], color = 'k', alpha = .5, label = "F_t_min");
        ax4[0, 2].axhline(self.horizonPlanner.F_t_lims[1], color = 'k', alpha = .5, label = "F_t_max");
        
        ax4[1, 2].axhline(np.rad2deg(self.horizonPlanner.alpha_lims[0]), color = 'k', alpha = .5, label = "alpha_min");
        ax4[1, 2].axhline(np.rad2deg(self.horizonPlanner.alpha_lims[1]), color = 'k', alpha = .5, label = "alpha_max");
        
        ax4[2, 2].axhline(np.rad2deg(self.horizonPlanner.beta_lims[0]), color = 'k', alpha = .5, label = "beta_min");
        ax4[2, 2].axhline(np.rad2deg(self.horizonPlanner.beta_lims[1]), color = 'k', alpha = .5, label = "beta_max");
        
        ax4[0, 2].set_ylim(0, 1.25*self.horizonPlanner.F_t_lims[1]);
        ax4[1, 2].set_ylim(1.25*np.rad2deg(self.horizonPlanner.alpha_lims[0]), 1.25*np.rad2deg(self.horizonPlanner.alpha_lims[1]));
        ax4[1, 2].set_ylim(1.25*np.rad2deg(self.horizonPlanner.beta_lims[0]), 1.25*np.rad2deg(self.horizonPlanner.beta_lims[1]));
        
        ax4[0, 2].legend(loc = 'best');
        ax4[1, 2].legend(loc = 'best');
        ax4[2, 2].legend(loc = 'best');
        