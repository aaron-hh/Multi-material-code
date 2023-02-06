This directory contains all the code added to perform cavitation modelling using Saurel model (two phase model).

eos_v2.cpp - This file contains the stiffened gas equation of state functions.
eos_v2.H - This file contains the function headers for eos_v2.cpp.
AmrLevelAdv.cpp - This originated 1-D advection code from Dr. Stephen Millmore has been modified to work on computing 2-D cavitation model. Changes can be tracked by following code labelled with '2020H'.
AmrLevelAdv.H - This file contains the function headers for AmrLevelAdv.cpp.
hllc.cpp - This file contains the functions for MUSCL Hancock HLLC solver.
hllc.H - This file contains the function headers for hllc.cpp.
CNS_bcfill.cpp - This file contains the external Dirichlet function defining the vibrating velocity equation.