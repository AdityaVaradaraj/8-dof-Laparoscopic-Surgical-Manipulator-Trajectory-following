Packages to unzip into src folder of catkin workspace:
 1) manipulator
 2) manipulator_moveit_config

I have submitted these folders/packages within the zip file alongwith the ppt and this README.txt 

Our aim is to make the robot follow 3 main trajectories:
1) Straight line to enter the wound or cut.
2) Arc in one plane to reorient the tool inside the wound.
3) Circle in a perpendicular plane to show that the tool can reach conventionally unreachable areas inside the wound to perform operations.

Algorithm: (moveit_joint_publisher.py)

We divide each trajectory into no. of points/steps and use Inverse Velocity Kinematics (q_dot = J_inverse * X_dot) to calculate joint velocities.

Then, we use q_current = q_previous + q_dot_current * dt 
similar to HW4 for each point in trajectory.

The, we use Moveit commands to pass in these q_current (, i.e, joint angle) values to the robot. 

Steps to Launch files:  

1) On one terminal, run the following command : 
	roslaunch manipulator_moveit_config demo.launch

This should load up the Moveit Rviz environment with our manipulator in it.

2) In Displays window, under "Planned Path", 
   a) Uncheck the checkbox next to "Loop Animation"
   b) Under Links --> base_link, select/check the checkbox next to "Show Axes"
   c) Under Links --> inst_link_3, select/check the checkbox next to "Show Trail"
   d) Under Panels (on the toolbar at the top), select / check the checkbox next to "Motion Planning - Trajectory Slider".
      This allows us to track the trajectory that the joint between the previous link and the small end-effector follows.
      Since we couldn't find any option to select the tip of the tip of the end-effector instead, this is a good enough 
      approximation since the end-effector is small in size.   
      
3) Now, in another terminal, launch the python file: 
	rosrun manipulator moveit_joint_publisher.py

4) See the trajectory being followed traced in the MoveIt Rviz environment.
