# DFF-Manipulation
Dynamic Flex-and-Flip Manipulation of Deformable Linear Objects


## Overview
This repository presents the implementation of flex-and-flip
manipulation tehcnique, suitable for grasping thin, flexible linear
objects lying on a flat surface. For example, we will show how page turning can be be performed with our methods.
Thus, the repository will contain:
1. Code to model deformable linear objects with minimum bending energy curves and two point contacts.
2. Detail fabrication of a soft pneumatic robotic gripper with two fingers suited for the task.
3. And the code to control the UR10 arm on which the gripper can be mounted

Here we assume that an AprilTag is printed on the edge of the page strip to know its position. Of course, this position can be hardcoded.

The following figures show the process of our flex-and-flip manipulation and how it is applied for a page turning task.

![figure2](media/fig2)

![figure1](media/fig1 | width=100)

## Usage
Here we expalin different elements of our pipeline.

### Modeling and Manipulating a Linear Deformable Object
The folder `modeling` 


