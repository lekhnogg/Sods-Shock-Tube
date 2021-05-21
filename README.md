# Shock-Tube-Simulation-Sod-
This code implicitly models fluid dynamics via a hydrodynamic simulation. 
We consider a tube filled with fluid. 
The fluid is separated by a diaphragm placed in the middle of the tube. Each partition contains a different density and pressure. 
By treating the fluid as an ideal gas, and considering more, smaller, discrete partitions,
we can model the interaction between the chambers once the diaphragm is lifted via shockwaves, rarefaction waves, 
and contact discontinuities between the discrete partitions. Sod's shock tube is a common test for the accuracy of computational  fluid codes (Wikipedia). 
This classic test problem is efficient in that while it only requires the modelling of a one dimensional system, we can study three different wave phenomena that emerge, 
which can be used to test the accuracy of this simulation results.

This code uses a Godunov based code to solve Sod's Shock tube problem.
