**Dewarp**

Lidar scans from a moving robot warp the world geometry causing problems for scan matching and map construction.

This project demonstrates a method to fix warped scans. By knowing the robot's linear and angular velocities a scan similar to what would be seen if the robot were standing still can be synthesized.

**Hacking**

This can be built and compiled using Visual Studio Code.  I've only tried it with GCC 5.4.1 / linux and it probably needs modifications to work in other environments.

This requires Eigen to build. 

Place a copy of Eigen in  ${workspaceFolder}/eigen34.