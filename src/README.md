
# Code

The src directory contains the whole source code that implements the fieldline 
tracing, interpolation and intersection testing algorithms.


The following code implements a TLAS (Top-Level-Accelerated-Structure)
behavior, that can contain multiple BLAS (Bottom-Level-Accelerated-Structure)
objects:

 - tlas.cpp
 - tlas.hpp

The code that implements the bicubic spline interpolation method:

 - bicubic.hpp
 - bicubic.cpp

The code that contains the RKF45 solver method (not actually used in the flt
code but used for testing the solver):

 - rkf45.cpp
 - rkf45.hpp

And finally the code that glues everything together via objects:

 - flt.cpp
 - flt.hpp

The code itself should be self-descriptive, therefore the best way is to dive
into the code to understand the inner workings.
