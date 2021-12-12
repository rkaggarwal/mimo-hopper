# -*- coding: utf-8 -*-
"""
Created on Sat May 29 16:09:27 2021

@author: raggarwal
"""

from casadi import *
import numpy as np
import matplotlib.pyplot as plt


"""
This class solves a single horizon given a dynamics model.  It uses Casadi
to manage the variables before solving with IPOPT, and defines everything 
for the optimization such as decision variables, objective function, and relevant
constraints.  This class does not contain the dynamics equations directly: instead,
it maintains a reference to a dynamics model and calls that dynamics model's 
derivative function to tie its collocation points together.

Definition of state and control vectors    

state: s = [x, y, z, xdot, ydot, zdot, q0, q1, q2, q3, body_omega_x, body_omega_y, body_omega_z]
control: u = [F_thrust, alpha, beta]

where [q0, q1, q2, q3] is a unit quaternion expressing the attitude of the body
with respect to the global inertial frame

x, y, and z are in the inertial frame (X points to the sky)
omega_x, omega_y and omega_z are with respect to the body frame

F_thrust is in N
alpha is the body y-axis servo angle (servo 1), rad
beta is the body z-axis servo angle (servo 2), rad
"""



class HorizonPlanner:

    def __init__(self, dynamicsModel,
                 no_of_horizon_steps = 50,
                 horizon_step_length_dt = .5,
                 F_t_lims = [2*9.81, 4*9.81],
                 alpha_lims = [-np.deg2rad(30), np.deg2rad(30)],
                 beta_lims = [-np.deg2rad(30), np.deg2rad(30)]):
            
        self.dynamicsModel = dynamicsModel;
        self.N = no_of_horizon_steps;
        self.dt = horizon_step_length_dt;
        self.horizon_length_t = self.N*self.dt;
        
        self.F_t_lims = F_t_lims;
        self.alpha_lims = alpha_lims;
        self.beta_lims = beta_lims;
        
        
        
    def planHorizon(self, s_init, s_final):
        """
        Plan a horizon from state-to-state
        Input: s_init, s_final
        Output: a mega array of z_maximal_horizon + u_horizon (maximal state and control evolution over the horizon)
        """
        
        ## SETUP ###
        self.opti = Opti();     

        ### DECISION VARIABLES ###
        S = self.opti.variable(13, self.N+1); # state
        U = self.opti.variable(3, self.N) # control effort evolution
        

        ### OBJECTIVE FUNCTION ###
        # Minimize quadratic cost on control effort
        cost = 0;
        Q = MX(diag([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])); # state weight bilinear form
        R = MX(diag([1, 1, 1])); # control weight bilinear form  
        
        Diag3 = MX(diag([1, 1, 1]));
        Diag4 = MX(diag([10, 10, 10, 10]));
        
        
        for k in range(0, self.N+1):
            
            cost += bilin(Q, S[:, k]-s_final, S[:, k]-s_final);
            
            if(k < self.N):
                cost += bilin(R, U[:, k], U[:, k]);

        self.opti.minimize(cost);  
        
        
        
        ### GENERAL CONSTRAINTS ###     
        
        self.opti.subject_to(S[0, :] >= 0); # X has to be above the ground at all times
        self.opti.subject_to(S[10, :] == 0); # no roll dynamics
        # self.opti.subject_to(DT >= 0);

       
        for k in range(self.N):
            
            # Euler Integration
            k1 = self.dynamicsModel.ds_dt(S[:, k], U[:, k]);
            s_next = S[:, k] + self.dt*k1
            
            
            # Re-normalize the attitude quaternion to have its 2-norm == 1 
            # (see http://ancs.eng.buffalo.edu/pdf/ancs_papers/2013/geom_int.pdf)
            
            s_next_q_norm = sqrt(s_next[6]**2 + s_next[7]**2 + s_next[8]**2 + s_next[9]**2);
            s_next[6] *= 1/s_next_q_norm;
            s_next[7] *= 1/s_next_q_norm;
            s_next[8] *= 1/s_next_q_norm;
            s_next[9] *= 1/s_next_q_norm;
            
            
            self.opti.subject_to(S[:, k+1] == s_next)
            

        ### CONTROL AND OTHER STATE CONSTRAINTS ###   
        
        # control effort constraints
        self.opti.subject_to(self.opti.bounded(self.F_t_lims[0],   U[0, :], self.F_t_lims[1]));
        self.opti.subject_to(self.opti.bounded(self.alpha_lims[0], U[1, :], self.alpha_lims[1]));
        self.opti.subject_to(self.opti.bounded(self.beta_lims[0],  U[2, :], self.beta_lims[1]));

        
        
        ### BOUNDARY CONDITIONS ###
     
        # Force the beginning of the horizon to coincide with the current
        # vehicle state as input into this overall function
        
        self.opti.subject_to(S[:, 0]  == s_init);
        self.opti.subject_to(S[:, -1]  == s_final);
        #self.opti.subject_to(S[3:13, -1] == s_final[3:13]); # and all the vel's & attitude the same
        # leave y and z free
        
        
        
        
        ### SOLVER INITIALIZATION ###
        
        # start off with just a linearly-interpolated kinematic path assumption
        # between the start and end states
        
        self.opti.set_initial(S[0, :], np.linspace(s_init[0], s_final[0], self.N+1));
        self.opti.set_initial(S[1, :], np.linspace(s_init[1], s_final[1], self.N+1));
        self.opti.set_initial(S[2, :], np.linspace(s_init[2], s_final[2], self.N+1));
        
        self.opti.set_initial(S[3, :], np.linspace(s_init[3], s_final[3], self.N+1));
        self.opti.set_initial(S[4, :], np.linspace(s_init[4], s_final[4], self.N+1));
        self.opti.set_initial(S[5, :], np.linspace(s_init[5], s_final[5], self.N+1));
        
        self.opti.set_initial(S[6, :], np.linspace(s_init[6], s_final[6], self.N+1));
        self.opti.set_initial(S[7, :], np.linspace(s_init[7], s_final[7], self.N+1));
        self.opti.set_initial(S[8, :], np.linspace(s_init[8], s_final[8], self.N+1));
        self.opti.set_initial(S[9, :], np.linspace(s_init[9], s_final[9], self.N+1));
        
        self.opti.set_initial(S[10, :], np.linspace(s_init[10], s_final[10], self.N+1));
        self.opti.set_initial(S[11, :], np.linspace(s_init[11], s_final[11], self.N+1));
        self.opti.set_initial(S[12, :], np.linspace(s_init[12], s_final[12], self.N+1));
                              


        # or if you have a kinematic reference, use that
        # self.opti.set_initial(S[0, :], xyz_ref[0, :]);
        # self.opti.set_initial(S[1, :], xyz_ref[1, :]);
        # self.opti.set_initial(S[2, :], xyz_ref[2, :]);
        
        # self.opti.set_initial(S[3, :], np.linspace(s_init[3], s_final[3], self.N+1));
        # self.opti.set_initial(S[4, :], np.linspace(s_init[4], s_final[4], self.N+1));
        # self.opti.set_initial(S[5, :], np.linspace(s_init[5], s_final[5], self.N+1));
        
        # self.opti.set_initial(S[6, :], q_flip_ref[0, :]);
        # self.opti.set_initial(S[7, :], q_flip_ref[1, :]);
        # self.opti.set_initial(S[8, :], q_flip_ref[2, :]);
        # self.opti.set_initial(S[9, :], q_flip_ref[3, :]);
        
        # self.opti.set_initial(S[10, :], np.linspace(s_init[10], s_final[10], self.N+1));
        # self.opti.set_initial(S[11, :], np.linspace(s_init[11], s_final[11], self.N+1));
        # self.opti.set_initial(S[12, :], np.linspace(s_init[12], s_final[12], self.N+1));
    

        #######################################################################
    
    
        ### SOLVE ###
    
    
    
        ### Solver Setup and Processing ###
        p_opts = {"expand": True}
        s_opts = {"print_level" : 5, "max_iter" : 1000};
        self.opti.solver("ipopt", p_opts, s_opts)
        sol = self.opti.solve()
    

    
        ### Process Results ###
        
        S_sol = sol.value(S)
        U_sol = sol.value(U)
        
        PHI_sol = np.zeros((self.N+1)); # roll angle
        THETA_sol = np.zeros((self.N+1)) # pitch angle
        PSI_sol = np.zeros((self.N+1)); # yaw angle
        
                 
        for i in range(0, self.N+1):
            # insert arithmetic here to convert from quaternion to roll/pitch/yaw
            q0 = S_sol[6, i];
            q1 = S_sol[7, i];
            q2 = S_sol[8, i];
            q3 = S_sol[9, i];
            
            PHI_sol[i] = np.arctan2(2*(q0*q1 + q2*q3), 1-2*(q1**2 + q2**2));
            THETA_sol[i] = np.arcsin(2*(q0*q2 - q3*q1));
            PSI_sol[i] = np.arctan2(2*(q0*q3 + q1*q2), 1-2*(q2**2 + q3**2));


        # Before appending our control efforts to the HORIZON_sol array, we need
        # to artificially increase their length by 1 for dimension matching
        # We do this at the end as just z.o.h. the last command (this doesn't 
        # affect the solution/optimization, it's just for diagnostics and convenience)
        
        
        U_sol = np.insert(U_sol, self.N-1, U_sol[:, -1], axis = 1);
        
   
        HORIZON_sol = np.vstack((S_sol, PHI_sol, THETA_sol, PSI_sol, U_sol));
        
        return HORIZON_sol; 
    