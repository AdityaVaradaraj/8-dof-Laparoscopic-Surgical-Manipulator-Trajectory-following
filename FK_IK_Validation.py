from sympy import *
import numpy as np
import math
from pprint import *
import matplotlib.pyplot as plt

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

#--------------- Inverse Kinematics for line -----------------------
N = 50          # no. of steps
t = 4            # time to complete
dt = t/N          # Time step
q_t = zeros(8,N+1)
q_t[:,0] = Matrix([[0],[0.253], [0.175],[0],[-1.744],[0],[0],[0]])
q_dot_t = zeros(8,N+1)

qa = Matrix([[0],[0.253], [0.175],[0],[-1.744],[0],[0],[0]])  
ef_pos = T_0_i[7].subs({"q1": qa[0], "q2": qa[1], "q3": qa[2], "q4": qa[3], "q5": qa[4], "q6": qa[5], "q7": qa[6],"q8": qa[7]})
# pprint(ef_pos)
home_pos = ef_pos[0:3, 3]
x_t0 = home_pos[0]
z_t0 = home_pos[2]
# --------------------------- Surgical instrument x,y,z for position inside body--------------------------------
qb = Matrix([[0],[0.253], [0.085],[0],[-1.744],[0],[0],[0]])  
ef_in_body_pos = T_0_i[7].subs({"q1": qb[0], "q2": qb[1], "q3": qb[2], "q4": qb[3], "q5": qb[4], "q6": qb[5], "q7": qb[6],"q8": qb[7]})
# pprint(ef_in_body_pos)
body_pos = ef_in_body_pos[0:3, 3]

x_dot = (body_pos[0]-home_pos[0])/4 # dx/dt
y_dot = 0 # No shift in y direction for line trajectory
z_dot = (body_pos[2]-home_pos[2])/4 # dz/dt

for i in range(N):
    X_dot = Matrix([[x_dot], [0], [z_dot], [0], [0], [0]])
    J_t = J.subs({"q1": q_t[0,i], "q2": q_t[1,i], "q3": q_t[2,i], "q4": q_t[3,i], "q5": q_t[4,i], "q6": q_t[5,i],"q7": q_t[6,i],"q8": q_t[7,i]})            
    J_inv = J_t.pinv()
        
    q_dot_t[:,i] = J_inv*X_dot
    q_t[:,i+1] = q_t[:,i] + q_dot_t[:,i]*dt


# -------------- Forward Kinematics ------------

X_dot_t = zeros(6, N+1)
X_t = zeros(6, N+1)
X_t[:, 0] = Matrix([[x_t0],[0],[z_t0],[0],[0],[0]])
#plt.axis([-0.150, 0.150, 0.530, 0.830])
plt.title("End-effector position using Velocity Kinematics")
plt.xlabel("X (in mts.)")
plt.ylabel("Z (in mts.)")
plt.scatter(x_t0, z_t0, color = 'blue')
plt.pause(0.05)
for i in range(N):
    J_t = J.subs({"q1": q_t[0,i], "q2": q_t[1,i], "q3" : q_t[2,i], "q4": q_t[3,i], "q5": q_t[4,i], "q6": q_t[5,i], "q7": q_t[6,i], "q8": q_t[7,i]})
    X_dot_t[:, i] = J_t*q_dot_t[:,i]
    X_t[:, i+1] = X_t[:,i] + X_dot_t[:,i]*dt
    x_t = X_t[0,i+1]
    z_t = X_t[2, i+1]
    plt.scatter(x_t, z_t, color = 'blue')
    plt.pause(0.05)
plt.show()