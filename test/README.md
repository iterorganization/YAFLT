
# Tests

The tests ensure the code is working as intended, even though that in some 
tests there is no actual check to perform, but are meant to see if the program
finishes with no errors

## Interpolation tests

In interpolation directory, the tests check the accuracy of the implementation
of the bicubic interpolation 2D algorithm, as well as it's performance and it's 
scalability. Actually the serial performance is not checked but is run to see
if the program finishes.

 - accuracy.cpp
   
    Tests the accuracy of the interpolation by intepolating an analytical 
    function (x - a)**2 * (y - b)**2 and seeing how many points differ by more 
    than .01%. It outputs a huge amount of text which can then be used to debug 
    in which regions the algorithm has issues (if you plot you should see an L
    band). Use the plot_points.py on the output of the program.

 - serial_performance.cpp
 
    Measures the time the algorithm computes N values in serial mode.

 - omp_performance.cpp
 
    Measures the time the algorithm computes N values in serial and parallel 
    mode. For evaluating performance scalability.

 - finite_difference_check.cpp 
 
    Checking the derivative tables if they are written exactly the same. Could 
    be used to see different schemas to be used to calculate derivative values.

## RKF45 tests

Tests the stability of the method by following a magnetic surface for a large 
number of steps. At each step evaluate if the new position is still on the 
magnetic surface or how much it differs from the initial flux value (should be 
constant). Of course the resolution of the input data affects the accuracy, 
nonetheless, the method itself should have some control over it.

Additionally there are two tests that contain the implementation of theRKF 4(5)
method and seeing if it is suitable for tracing fieldlines.

 - test_stability.cpp
 - test_rkf45_method.cpp
 - test_rkf45_method_stability.cpp

## Creation and deletion of objects

This checks for memory leaks and creation/deletion of objects (for finding if 
a segfault is lurking somewhere).

 - create_and_delete_objects.cpp

## Embree tests

This tests loading mesh data to Embree and checking if intersection tests work 
on a simple case.

 - test_intersection.cpp
 - load_mesh.cpp