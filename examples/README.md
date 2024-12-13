
# Examples

This directory contains examples of how to use the code to obtain fieldlines or
how to run a simple FLT case.

The examples also output VTK ASCII files that can be visualized in ParaView.

The magnetic data is stored in the header files.

flt1.cpp obtains points of a fieldline and outputs them in a VTK file. 

flt2.cpp obtains a fieldline by followining it from the barycenter of a 
triangle and outputing the points in a VTK file. 

flt3.cpp is same as flt2, except that in this case, based on the input normals we
correctly determine the direction we should follow the fieldlines.

And inres1_example.cpp runs the inres1 case. This is a limiter equilibrium, 
touching FWP #4. The input data is included in headers:
 - data_for_shadow.h
   Shadowing geometry data, used for intersection tests.
 - data_for_target.h
   Target geometry data.
 - data_for_spline.h
   Equilibrium data, aka 'EQ3.eqdsk'

The output is a vtk file, containing the connection lengths information and 
the intersection mask arrays.