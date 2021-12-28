from sympy import *
import numpy as np
import math
from pprint import *
from mpl_toolkits.mplot3d import Axes3D
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

# print("Final Transformation Matrix: ")
# pprint(T_0_i[7])
# print()

x = []
y = []
z = []
# ------------------------ Workspace Study -------------------------
for qw1 in np.arange(-185*np.pi/180, 185*np.pi/180 + np.pi/4, np.pi/4):
    for qw2 in np.arange(50*np.pi/180, -np.pi/2 - 35*np.pi/180, -35*np.pi/180):
        for qw3 in np.arange(60*np.pi/180, -138*np.pi/180 - 49.5*np.pi/180, -49.5*np.pi/180):
            qw_values = {"q1": qw1, "q2": qw2, "q3": qw3, "q4": 0, "q5": 0, "q6": 0, "q7": 0, "q8": 0}
            T = T_0_i[7].subs(qw_values)
            x.append(T[0, 3])
            y.append(T[1, 3])
            z.append(T[2,3])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(z, dtype=float), marker='o')

ax.set_xlabel('X (mts.)')
ax.set_ylabel('Y (mts.)')
ax.set_zlabel('Z (mts.)')

plt.show()