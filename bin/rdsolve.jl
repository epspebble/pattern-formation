"""
Solution a generic reaction-diffusion (RD) system using Method of Lines (MoL), with the DifferentialEquations
package on Julia.

This program should define:

1) a struct that represents a RD system in its mathematical entirety, should contain PDE parameters, 
and at least initial and boundary conditions, and then also full temporal-spatial simulation data, 
when a time-stepper is called to act on this system. 

2) a small time-stepper function that runs the simulation at a user-supplied amount of time change, 
using the DifferentialEquations Julia module, by simply producing an ODEProblem object using a method of line 
approach as an intermediate step and calls solve(), with the default choice of CVODE_BDF(linear_solver:GMRES)) 
as the ODE solver.

3) a load/save routine for transferring the state of the RD system to and from a netCDF file on the disk.

4) basic visualization functions

"""
