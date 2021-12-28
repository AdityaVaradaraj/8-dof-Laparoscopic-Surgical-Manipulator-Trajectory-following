#!/usr/bin/env python
from __future__ import print_function
from six.moves import input

import sys
import copy
import rospy
from sympy import *
import numpy as np
import math
from pprint import *
import matplotlib.pyplot as plt
import moveit_commander
import moveit_msgs.msg
import geometry_msgs.msg
def cross_product(Z, O):
    C = zeros(3, 1)
    C[0] = Z[1]*O[2] - Z[2]*O[1]
    C[1] = -(Z[0]*O[2] - Z[2]*O[0])
    C[2] = Z[0]*O[1] - Z[1]*O[0]
    return C


q1, q2, q3, q4, q5, q6, q7, q8 = symbols('q1, q2, q3, q4, q5, q6, q7, q8')

# ------------------------- ALL DISTANCES IN MTS. ----------------------------
# ----------------------------- DH PARAMETERS --------------------------------
a_i = [0.330, 0.1150, 0, 0, 0, 0, 0.058, 0.020]
d_i = [0.645, 0, 0, 1.220,  0, 0.437,  0, 0]
alpha_i = [np.pi/2, 0, np.pi/2,-np.pi/2, np.pi/2, np.pi/2, np.pi/2, 0]
theta_i = [np.pi + q1, np.pi/2+q2, q3, q4, q5, np.pi/2+q6, np.pi/2+q7, q8] 

# dof --> degrees of freedom
dof = 8

# --------------------------- Transformation Matrices -----------------------
T_0_1 = eye(4)
T_0_2 = eye(4)
T_0_3 = eye(4)
T_0_4 = eye(4)
T_0_5 = eye(4)
T_0_6 = eye(4)
T_0_7 = eye(4)
T_0_n = eye(4)

T_0_i = [T_0_1, T_0_2, T_0_3, T_0_4, T_0_5, T_0_6, T_0_7, T_0_n] 
Z_i = [zeros(3,1), zeros(3,1), zeros(3,1), zeros(3,1), zeros(3,1), zeros(3,1), zeros(3,1), zeros(3, 1)]
O_i = [zeros(3,1), zeros(3,1), zeros(3,1), zeros(3,1), zeros(3,1), zeros(3,1), zeros(3,1), zeros(3, 1), zeros(3,1)]
Z_i[0] = Matrix([[0],[0],[1]])
J = zeros(6,8)


for j in range(dof):
    A_i = Matrix([[cos(theta_i[j]), -sin(theta_i[j])*cos(alpha_i[j]), sin(theta_i[j])*sin(alpha_i[j]), a_i[j]*cos(theta_i[j])], 
                [sin(theta_i[j]), cos(theta_i[j])*cos(alpha_i[j]), -cos(theta_i[j])*sin(alpha_i[j]), a_i[j]*sin(theta_i[j])],
                [0, sin(alpha_i[j]), cos(alpha_i[j]), d_i[j]],
                [0, 0, 0, 1]])
    


    if(j == 0):
        T_0_i[j] = A_i
    else:
        T_0_i[j] = T_0_i[j-1]*A_i

    if j != dof-1:
        Z_i[j+1] = T_0_i[j][0:3, 2]
    O_i[j+1] = T_0_i[j][0:3, 3]

# ---------------------------   Jacobian Calculation --------------------------------
    
for j in range(dof):
    J[0:3, j] = cross_product(Z_i[j], O_i[dof] - O_i[j])
    J[3:, j] = Z_i[j]

# Joint Angles for home and inside body was obtained by playing around with angles in Gazebo hospital_scene.world

# -------------------------- Trajectory 1: Entering the body using line trajectory ----------------------

# --------------------------- Surgical instrument x,y,z for home position--------------------------------
q = Matrix([[0],[0.253], [0.175],[0],[-1.744],[0],[0],[0]])  
ef_pos = T_0_i[7].subs({"q1": q[0], "q2": q[1], "q3": q[2], "q4": q[3], "q5": q[4], "q6": q[5], "q7": q[6],"q8": q[7]})
# pprint(ef_pos)
home_pos = ef_pos[0:3, 3]

# --------------------------- Surgical instrument x,y,z for position inside body--------------------------------
q1 = Matrix([[0],[0.253], [0.085],[0],[-1.744],[0],[0],[0]])  
ef_in_body_pos = T_0_i[7].subs({"q1": q1[0], "q2": q1[1], "q3": q1[2], "q4": q1[3], "q5": q1[4], "q6": q1[5], "q7": q1[6],"q8": q1[7]})
# pprint(ef_in_body_pos)
body_pos = ef_in_body_pos[0:3, 3]

x_dot = (body_pos[0]-home_pos[0])/4 # dx/dt
y_dot = 0 # No shift in y direction for line trajectory
z_dot = (body_pos[2]-home_pos[2])/4 # dz/dt



try:
    from math import pi, tau, dist, fabs, cos
except:  # For Python 2 compatibility
    from math import pi, fabs, cos, sqrt

    tau = 2.0 * pi

    def dist(p, q):
        return sqrt(sum((p_i - q_i) ** 2.0 for p_i, q_i in zip(p, q)))

# ---------------- MoveIt Initializations --------------------------
from std_msgs.msg import String
from moveit_commander.conversions import pose_to_list

moveit_commander.roscpp_initialize(sys.argv)
rospy.init_node("manipulator_joint_publisher", anonymous=True)

robot = moveit_commander.RobotCommander()

scene = moveit_commander.PlanningSceneInterface()

group_name = "manipulator"
move_group = moveit_commander.MoveGroupCommander(group_name)



# ---------------------------   Inverse Velocity Kinematics for Line Trajectory 1 --------------------------------

N = 4          # no. of steps
t = 4            # time to complete
dt = t/N          # Time step
q_t = zeros(8,N+1)
q_t[:,0] = Matrix([[0],[0.253], [0.175],[0],[-1.744],[0],[0],[0]]) # Initial Start configuration 
# Obtained by playing around with angles in Gazebo using rosservice call /manipulator/SetModelConfiguration

q_dot_t = zeros(8,N+1)
joint_goal = move_group.get_current_joint_values()
joint_goal[0] = float(q_t[0, 0] - 21*pi/180)
joint_goal[1] = float(q_t[1, 0] - 30*pi/180)
joint_goal[2] = float(q_t[2, 0] - 3*pi/180)
joint_goal[3] = float(q_t[3, 0] + 24*pi/180)
joint_goal[4] = float(q_t[4, 0] + pi/2 - 62*pi/180)
joint_goal[5] = float(q_t[5, 0] + 0.0)  
joint_goal[6] = float(q_t[6, 0] + 107*pi/180)
joint_goal[7] = float(q_t[7, 0] - 108*pi/180)
move_group.go(joint_goal, wait=True)

move_group.stop()


for i in range(N):
    X_dot = Matrix([[x_dot], [0], [z_dot], [0], [0], [0]])
    J_t = J.subs({"q1": q_t[0,i], "q2": q_t[1,i], "q3": q_t[2,i], "q4": q_t[3,i], "q5": q_t[4,i], "q6": q_t[5,i],"q7": q_t[6,i],"q8": q_t[7,i]})            
    J_inv = J_t.pinv()
        
    q_dot_t[:,i] = J_inv*X_dot
    q_t[:,i+1] = q_t[:,i] + q_dot_t[:,i]*dt
    joint_goal = move_group.get_current_joint_values()
    
    # We add the error to the q values before passing to MoveIt to convert from our DH space to MoveIt space
    # E.g.: -21 deg is error compensation for difference between 
    # angle values acc. to axes taken for DH and angles to corresponding configuration in MoveIt
    
    joint_goal[0] = float(q_t[0, i+1] - 21*pi/180) 
    joint_goal[1] = float(q_t[1, i+1] - 30*pi/180)
    joint_goal[2] = float(q_t[2, i+1] - 3*pi/180)
    joint_goal[3] = float(q_t[3, i+1] + 24*pi/180)
    joint_goal[4] = float(q_t[4, i+1] + pi/2 - 62*pi/180)
    joint_goal[5] = float(q_t[5, i+1] + 0.0)  
    joint_goal[6] = float(q_t[6, i+1] + 107*pi/180)
    joint_goal[7] = float(q_t[7, i+1] - 108*pi/180)
    move_group.go(joint_goal, wait=True)
move_group.stop()    

# ---------------------------------- Reorienting tool inside the wound ----------------------------------------
# --------------------------- Inverse Velocity Kinematics for Arc Trajectory 2 --------------------------------

N1 = 6          # no. of steps
t = 2            # time to complete
dt = t/N1          # Time step
q_t_1 = zeros(8,N1+1)
th_dot = np.pi/(2*t)
q_dot_t_1 = zeros(8,N1+1)
q_t_1[:,0] = q_t[:, N] 
r1 = 0.078 # Radius of Arc in Mts.
th = 0
dth = np.pi/(2*N1)
for i in range(N1):
    X_dot = Matrix([[-th_dot*r1*cos(th)], [0], [th_dot*r1*sin(th)], [0], [0], [0]])
    J_t = J.subs({"q1": q_t_1[0,i], "q2": q_t_1[1,i], "q3": q_t_1[2,i], "q4": q_t_1[3,i], "q5": q_t_1[4,i], "q6": q_t_1[5,i],"q7": q_t_1[6,i],"q8": q_t_1[7,i]})            
    J_inv = J_t.pinv()
        
    q_dot_t_1[:,i] = J_inv*X_dot
    q_t_1[:,i+1] = q_t_1[:,i] + q_dot_t_1[:,i]*dt
    joint_goal = move_group.get_current_joint_values()

    # We add the error to the q values before passing to MoveIt to convert from our DH space to MoveIt space
    # E.g.: -21 deg is error compensation for difference between 
    # angle values acc. to axes taken for DH and angles to corresponding configuration in MoveIt

    joint_goal[0] = float(q_t_1[0, i+1] - 21*pi/180)
    joint_goal[1] = float(q_t_1[1, i+1] - 30*pi/180)
    joint_goal[2] = float(q_t_1[2, i+1] - 3*pi/180)
    joint_goal[3] = float(q_t_1[3, i+1] + 24*pi/180)
    joint_goal[4] = float(q_t_1[4, i+1] + pi/2 - 62*pi/180)
    joint_goal[5] = float(q_t_1[5, i+1] + 0.0)  
    joint_goal[6] = float(q_t_1[6, i+1] + 107*pi/180)
    joint_goal[7] = float(q_t_1[7, i+1] - 108*pi/180)
    move_group.go(joint_goal, wait=True)
    
    th += dth
move_group.stop()


# ----------------------- Reaching conventionally inaccessible areas inside the wound ----------------------------
# --------------------------- Inverse Velocity Kinematics for Circle Trajectory 2 --------------------------------

N2 = 15          # no. of steps
t = 4            # time to complete
dt = t/N1          # Time step
q_t_2 = zeros(8,N2+1)
th_dot = 2*np.pi/(t)
q_dot_t_2 = zeros(8,N2+1)
q_t_2[:,0] = q_t_1[:, N1] 
r2 = 0.078 # Radius of Circle in Mts.
th = 0
dth = 2*np.pi/(N2)
for i in range(N2):
    X_dot = Matrix([[0], [-th_dot*r1*cos(th)], [th_dot*r1*sin(th)], [0], [0], [0]])
    J_t = J.subs({"q1": q_t_2[0,i], "q2": q_t_2[1,i], "q3": q_t_2[2,i], "q4": q_t_2[3,i], "q5": q_t_2[4,i], "q6": q_t_2[5,i],"q7": q_t_2[6,i],"q8": q_t_2[7,i]})            
    J_inv = J_t.pinv()
        
    q_dot_t_2[:,i] = J_inv*X_dot
    q_t_2[:,i+1] = q_t_2[:,i] + q_dot_t_2[:,i]*dt
    joint_goal = move_group.get_current_joint_values()

    # We add the error to the q values before passing to MoveIt to convert from our DH space to MoveIt space
    # E.g.: -21 deg is error compensation for difference between 
    # angle values acc. to axes taken for DH and angles to corresponding configuration in MoveIt

    joint_goal[0] = float(q_t_2[0, i+1] - 21*pi/180)
    joint_goal[1] = float(q_t_2[1, i+1] - 30*pi/180)
    joint_goal[2] = float(q_t_2[2, i+1] - 3*pi/180)
    joint_goal[3] = float(q_t_2[3, i+1] + 24*pi/180)
    joint_goal[4] = float(q_t_2[4, i+1] + pi/2 - 62*pi/180)
    joint_goal[5] = float(q_t_2[5, i+1] + 0.0)  
    joint_goal[6] = float(q_t_2[6, i+1] + 107*pi/180)
    joint_goal[7] = float(q_t_2[7, i+1] - 108*pi/180)
    move_group.go(joint_goal, wait=True)
    
    th += dth
move_group.stop()