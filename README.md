countour2D_CUDA
===============

This is code for program, which calculates 2d contour points of polygonal data in CUDA C
Contour calculation of 2D polygonal data using CUDA

Recently I had a project for contour calculation, and tried to make it on NVIDIA graphic card, because it is one of the fastest platform nowadays on ordinary machines.
There are many computational geometry algorithms for making contour calculation of 2D point data, but almost none for making such calculations for array of polygons.
Provided algorithm uses brute force for obtaining result but works ok when used of GPU because of many multithreaded calculation.
Algorithm is very simple:
<ul>
<li>Detection of all figures' lines intersections</li>
<li>Calculation of all points, which lay on the edge between square "inside" contour and rest space, which is made with even more brute force: around each point program generates circle with Epsilon radius and with AngleStep step, and each point of generated circle is checked if it lays on even one of polygons of scene. If circle contains sector, which does not belong to any polygon, then it means, that this is edge point and, of course, point of contour.</li>
</ul>
Program works well, but can be optimized.
Hope you will be able to use, because I made it as simple as possible to compile.
