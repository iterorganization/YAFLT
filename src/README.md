
# Code

The src directory contains the whole source code that implements the fieldline 
tracing, interpolation and intersection testing algorithms.


The code that contains that makes a thin specialized layer for intersection 
tests on a triangular mesh via Embree is as following:

 - accell_embree.cpp
 - accell_embree.hpp

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
